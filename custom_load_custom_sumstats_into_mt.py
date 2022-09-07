#!/usr/bin/env python3

__author__ = 'rahulg modifying konradk'

import sys
import argparse
from pprint import pprint
from datetime import date
from collections import Counter
from tqdm import tqdm
from ukbb_pan_ancestry import *
from ukbb_pan_ancestry.load_all_results import apply_qc
from ukbb_pan_ancestry.saige_pan_ancestry_custom import *


def get_all_valid_variant_results_ht_paths(pop, suffix):
    results_dir = f'{root}/result/{suffix}/{pop}'
    all_phenos_dir = hl.hadoop_ls(results_dir)
    all_variant_outputs = get_files_in_parent_directory(all_phenos_dir)
    return all_variant_outputs


def get_pheno_dict(suffix):
    pheno_ht = hl.read_matrix_table(get_custom_ukb_pheno_mt_path(suffix)).cols()
    pheno_dict = create_broadcast_dict(pheno_ht.key)
    return pheno_dict


def custom_get_variant_results_path(pop: str, suffix: str, extension: str = 'mt'):
    return f'{root}/sumstats/{suffix}/mt/results_{pop}.{extension}'


def custom_unify_saige_ht_schema(ht, patch_case_control_count: str = ''):
    """

    :param Table ht:
    :param str patch_case_control_count: Path to file (hack to get cases and controls back if loading later)
    :return:
    :rtype: Table
    """
    assert ht.head(1).annotation.collect()[0] is None, f'failed at {patch_case_control_count}'
    if 'AF.Cases' not in list(ht.row):
        ht = ht.select('AC_Allele2', 'AF_Allele2', 'MissingRate', 'N', 'BETA', 'SE',
                       **{'p.value.NA': hl.null(hl.tfloat64), 'Is.SPA.converge': hl.null(hl.tint32),
                          'var': ht.var, 'AF.Cases': hl.null(hl.tfloat64),
                          'AF.Controls': hl.null(hl.tfloat64), 'Pvalue': ht.Pvalue,
                          'gene': hl.or_else(ht.gene, ''), 'annotation': hl.or_else(ht.annotation, '')})
    else:
        ht = ht.select('AC_Allele2', 'AF_Allele2', 'MissingRate', 'N', 'BETA', 'SE',
                       'p.value.NA', 'Is.SPA.converge', 'var', 'AF.Cases',
                       'AF.Controls', 'Pvalue', gene=hl.or_else(ht.gene, ''), annotation=hl.or_else(ht.annotation, ''))
    
    ht = ht.annotate(BETA = hl.float64(ht.BETA), SE = hl.float64(ht.SE))

    # ht2 = ht.head(1)
    # pheno_key_dict = dict(ht2.aggregate(hl.agg.take(ht2.key, 1)[0]))
    # if patch_case_control_count:
    #     if not ht.n_cases.collect()[0]:
    #         directory, tpc, _ = patch_case_control_count.rsplit('/', 2)

    #         pheno_results_dir = get_pheno_output_path(directory, pheno_key_dict, '', legacy=False)
    #         prefix = get_results_prefix(pheno_results_dir, pheno_key_dict, '{chrom}', 1, legacy=False)
    #         saige_log = f'{prefix}.variant.log'
    #         cases, controls = custom_get_cases_and_controls_from_log(saige_log)
    #         print(f'Patched pheno: {tpc}. Got {cases} cases and {controls} controls.')
    #         if cases == -1: cases = hl.null(hl.tint)
    #         if controls == -1: controls = hl.null(hl.tint)
    #         ht = ht.annotate_globals(n_cases=cases, n_controls=controls)
    if 'heritability' not in list(ht.globals):
        ht = ht.annotate_globals(heritability=hl.null(hl.tfloat64))
    if 'saige_version' not in list(ht.globals):
        ht = ht.annotate_globals(saige_version=hl.null(hl.tstr))
    return ht


def custom_patch_mt_keys(mt):
    mt = mt.key_cols_by(**{x: hl.case(missing_false=True)
                        .default(mt[x])
                           for x in PHENO_KEY_FIELDS})
    return mt


def check_and_annotate_with_dict(mt, input_dict, dict_key_from_mt, axis='cols'):
    direction = mt.col if axis == 'cols' else mt.row
    annotation = {}
    for key, value in list(input_dict.values()[0].items()):
        if key in list(direction):
            annotation[key] = hl.case().when(
                hl.is_defined(mt[key]), mt[key]).when(
                hl.is_defined(input_dict.get(dict_key_from_mt)), input_dict.get(dict_key_from_mt)[key]).default(
                hl.null(value.dtype)
            )
        else:
            annotation[key] = hl.if_else(hl.is_defined(input_dict.get(dict_key_from_mt)),
                                         input_dict.get(dict_key_from_mt)[key],
                                         hl.null(value.dtype))
    if axis == 'cols':
        mt = mt.annotate_cols(**annotation)
    else:
        mt = mt.annotate_rows(**annotation)
    return mt


def generate_sumstats_mt(all_variant_outputs, pheno_dict, temp_dir, inner_mode = '_read_if_exists',
                         checkpoint: bool = False):
    row_keys = ['locus', 'alleles', 'gene', 'annotation']
    col_keys = PHENO_KEY_FIELDS

    all_hts = [custom_unify_saige_ht_schema(hl.read_table(x), patch_case_control_count=x) for x in tqdm(all_variant_outputs)]
    mt = join_pheno_hts_to_mt(all_hts, row_keys, col_keys, temp_dir=temp_dir,
                              inner_mode=inner_mode, repartition_final=20000)

    if checkpoint:
        mt = mt.checkpoint(f'{temp_dir}/staging.mt', **{inner_mode: True})
    mt = custom_patch_mt_keys(mt)
    key = mt.col_key.annotate(phenocode=format_pheno_dir(mt.phenocode))
    mt = check_and_annotate_with_dict(mt, pheno_dict, key)
    if mt.inv_normalized.dtype == hl.tstr:
        mt = mt.annotate_cols(inv_normalized=hl.bool(mt.inv_normalized))
    key = mt.col_key.annotate(phenocode=format_pheno_dir(mt.phenocode))
    mt = mt.filter_cols(mt.phenocode != "").drop('Is.SPA.converge', 'AC_Allele2')
    mt = mt.key_rows_by('locus', 'alleles')
    return mt


def write_full_mt(suffix, temp_dir, overwrite):
    mts = []
    for pop in POPS:
        mt = hl.read_matrix_table(custom_get_variant_results_path(pop, suffix, 'mt')).annotate_cols(pop=pop)
        mt = mt.annotate_cols(_logged=hl.agg.any(mt.Pvalue < 0))
        mt = mt.annotate_entries(Pvalue=hl.if_else(mt._logged, mt.Pvalue, hl.log(mt.Pvalue))).drop('_logged')

        mt = apply_qc(mt)
        mt = custom_patch_mt_keys(mt)
        mt = reannotate_cols(mt, suffix)
        mt = re_colkey_mt(mt)
        mt = mt.select_cols(pheno_data=mt.col_value)
        mt = mt.select_entries(summary_stats=mt.entry)
        mts.append(mt)

    full_mt = mts[0]
    for mt in mts[1:]:
        full_mt = full_mt.union_cols(mt, row_join_type='outer')
    
    full_mt = full_mt.checkpoint(f'{temp_dir}/staging_full.mt', overwrite=True)
    full_mt = full_mt.collect_cols_by_key()
    full_mt = full_mt.annotate_cols(
        pheno_data=full_mt.pheno_data.map(lambda x: x.drop(*PHENO_COLUMN_FIELDS)),
        **{x: full_mt.pheno_data[x][0] for x in PHENO_DESCRIPTION_FIELDS},
        **{f'n_cases_full_cohort_{sex}': full_mt.pheno_data[f'n_cases_{sex}'][0]
           for sex in ('both_sexes', 'females', 'males')}
    )
    full_mt.write(custom_get_variant_results_path('full', suffix), overwrite)
    print('Pops per pheno:')
    pprint(dict(Counter(full_mt.aggregate_cols(hl.agg.counter(hl.len(full_mt.pheno_data))))))


def reannotate_cols(mt, suffix):
    pheno_dict = get_pheno_dict(suffix)
    key = get_modified_key(mt)
    mt = check_and_annotate_with_dict(mt, pheno_dict, key)
    return mt


def recode_single_pheno_struct_to_legacy_path(pheno_struct):
    trait_type = pheno_struct.trait_type
    phenocode = pheno_struct.phenocode
    if trait_type == 'icd10':
        trait_type = 'icd_all'
        coding = 'icd10'
    elif trait_type == 'phecode':
        coding = pheno_struct.pheno_sex
    elif trait_type == 'biomarkers':
        coding = pheno_struct.phenocode
    else:
        if phenocode == 'whr':
            coding = 'whr'
        else:
            coding = pheno_struct.coding if pheno_struct.coding else pheno_struct.modifier
    return f'{trait_type}-{format_pheno_dir(phenocode)}-{coding}'


def main(args):
    hl.init(default_reference='GRCh37', log='/combine_results.log', branching_factor=8)

    inner_mode = 'overwrite' if args.overwrite else '_read_if_exists'
    pops = args.pops.split(',') if args.pops else POPS

    if args.run_basic_load:
        for pop in pops:
            all_variant_outputs = get_all_valid_variant_results_ht_paths(pop, args.suffix)
            pheno_dict = get_pheno_dict(args.suffix)

            mt = generate_sumstats_mt(all_variant_outputs, pheno_dict, f'{temp_bucket}/{args.suffix}/{pop}/variant', inner_mode)
            mt.write(custom_get_variant_results_path(pop, args.suffix, 'mt'), overwrite=args.overwrite)

    if args.run_additional_load:
        today = date.today().strftime("%y%m%d")
        for pop in pops:
            # data will be drawn from --suffix2 and merged into --suffix
            all_variant_outputs = get_all_valid_variant_results_ht_paths(pop, args.suffix2)
            pheno_dict = get_pheno_dict(args.suffix2)

            all_keys = hl.read_matrix_table(custom_get_variant_results_path(pop, args.suffix, 'mt')).col_key.collect()
            loaded_phenos = {
                recode_single_pheno_struct_to_legacy_path(x) for x in all_keys
            }
            loaded_phenos |= {f'{x.trait_type}-{format_pheno_dir(x.phenocode)}-{x.pheno_sex}-{x.coding}-{x.modifier}' for x in all_keys}

            def _matches_any_pheno(pheno_path, phenos_to_match):
                return any(x for x in phenos_to_match if f'/{x}/variant_results.ht' in pheno_path)

            if args.force_reload:
                pheno_matches = set(args.force_reload.split(','))
                all_variant_outputs = [x for x in all_variant_outputs if _matches_any_pheno(x, pheno_matches)]
            else:
                all_variant_outputs = [x for x in all_variant_outputs if not _matches_any_pheno(x, loaded_phenos)]

                if args.load_only:
                    pheno_matches = set(args.load_only.split(','))
                    if '' in pheno_matches:
                        print('WARNING: Empty string in pheno_matches. Might reload more than expected')
                    all_variant_outputs = [x for x in all_variant_outputs if _matches_any_pheno(x, pheno_matches)]

            print(f'Loading {len(all_variant_outputs)} additional HTs...')
            if len(all_variant_outputs) < 20:
                print(all_variant_outputs)
            if args.dry_run:
                continue
            if not len(all_variant_outputs):
                continue

            mt = generate_sumstats_mt(all_variant_outputs, pheno_dict,
                                      f'{temp_bucket}/{args.suffix}/{pop}/variant_{today}', inner_mode, checkpoint=True)

            original_mt = hl.read_matrix_table(custom_get_variant_results_path(pop, args.suffix, 'mt'))
            original_mt = original_mt.checkpoint(f'{temp_bucket}/{args.suffix}/{pop}/variant_before_{today}.mt', overwrite=True)
            if args.force_reload:
                original = original_mt.count_cols()
                original_mt = original_mt.filter_cols(
                    hl.literal(pheno_matches).contains(hl.delimit(list(original_mt.col_key.values()), '-')), keep=False)
                print(f'\n\nGoing from {original} to {original_mt.count_cols()}...\n\n')
            mt = re_colkey_mt(mt.annotate_cols(pop=pop))
            original_mt = re_colkey_mt(original_mt.annotate_cols(pop=pop))
            mt = original_mt.union_cols(mt, row_join_type='outer')
            mt.write(custom_get_variant_results_path(pop, args.suffix, 'mt'), overwrite=args.overwrite)

    if args.run_combine_load:
        write_full_mt(args.suffix, f'{temp_bucket}/{args.suffix}/full', args.overwrite)


def re_colkey_mt(mt):
    mt = mt.key_cols_by().select_cols(*PHENO_KEY_FIELDS, *PHENO_COLUMN_FIELDS, *PHENO_GWAS_FIELDS, 'log_pvalue', 'pop')
    return mt.key_cols_by(*PHENO_KEY_FIELDS)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--suffix', help='Loads everything for a particular suffix.', required=True)
    parser.add_argument('--suffix2', help='Use with additional load. Will append items from this to MTs for suffix.')
    parser.add_argument('--run_basic_load', help='Overwrite everything', action='store_true')
    parser.add_argument('--run_additional_load', help='Overwrite everything', action='store_true')
    parser.add_argument('--run_combine_load', help='Overwrite everything', action='store_true')
    parser.add_argument('--dry_run', help='Overwrite everything', action='store_true')
    parser.add_argument('--load_only', help='Comma-separated list of trait_type-pheno-coding to run'
                                            '(e.g. continuous-50-irnt,icd_all-E10-icd10 )')
    parser.add_argument('--force_reload', help='Comma-separated list of trait_type-pheno-coding to force reload'
                                            '(e.g. continuous-50-irnt,icd_all-E10-icd10 )')
    parser.add_argument('--pops', help='comma-separated list')
    args = parser.parse_args()

    main(args)
