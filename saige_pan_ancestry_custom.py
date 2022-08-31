#!/usr/bin/env python3

import sys
import copy
import argparse
import logging
from datetime import date
logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s", level='INFO', filename='saige_pipeline.log')

from gnomad.utils import slack
from ukb_common import *
import time
import re

from ukbb_pan_ancestry import *
from ukb_common.utils.saige_pipeline import *

logger = logging.getLogger("saige_pan_ancestry_custom")
logger.addHandler(logging.StreamHandler(sys.stderr))
bucket_vcf = 'gs://ukb-diverse-pops'
root_vcf = f'{bucket_vcf}/results'
bucket = 'gs://rgupta-assoc'
root = f'{bucket}/saige_gwas'
pheno_folder = f'{bucket}/phenotype'

DESCRIPTION_PATH = f'{pheno_folder}/field_descriptions.tsv'

HAIL_DOCKER_IMAGE = 'gcr.io/1/hail_utils:6.1'
PHENO_DOCKER_IMAGE = 'gcr.io/ukbb-diversepops-neale/rgupta_hail_utils'
SAIGE_DOCKER_IMAGE = 'gcr.io/ukbb-diversepops-neale/saige:1.1.5'
QQ_DOCKER_IMAGE = 'konradjk/saige_qq:0.2'

SCRIPT_DIR = '/ukb_common/saige'
# TODO: add binary_trait annotation to input table and remove this:
saige_pheno_types = {
    'continuous': 'quantitative',
    'biomarkers': 'quantitative',
    'categorical': 'binary',
    'icd': 'binary',
    'icd_first_occurrence': 'binary',
    'icd_all': 'binary',
    'phecode': 'binary',
    'prescriptions': 'binary'
}


def get_custom_ukb_pheno_mt_path(suffix):
    return f'{pheno_folder}/mt/phenotype_{suffix}.mt'


def get_custom_phenotype_summary_backup_path(suffix, curdate):
    return f'{pheno_folder}/summary/all_pheno_summary_{suffix}_before_{curdate}.txt.bgz'


def get_custom_phenotype_summary_path(suffix, data_type: str, extension = 'ht'):
    return f'{pheno_folder}/summary/phenotype_{data_type}_{suffix}.{extension}'


def get_custom_munged_pheno_path(suffix):
    return f'{pheno_folder}/mt/munged/munged_raw_phenotype_{suffix}.mt'


def get_custom_ukb_pheno_mt(suffix, pop: str = 'all'):
    mt = hl.read_matrix_table(get_custom_ukb_pheno_mt_path(suffix))
    # mt = mt.annotate_rows(**get_ukb_meta(key_type=hl.tint32)[mt.row_key])
    mt = mt.annotate_rows(**get_covariates(key_type=hl.int32)[mt.row_key])
    if pop != 'all':
        mt = mt.filter_rows(mt.pop == pop)
    return mt


def custom_get_phenos_to_run(suffix, pop, limit, specific_phenos,
                             skip_case_count_filter, sex_stratified):
    ht = hl.read_table(get_custom_phenotype_summary_path(suffix))
    ht = ht.filter(ht.pop == pop)
    min_cases = MIN_CASES_EUR if pop == 'EUR' else MIN_CASES

    criteria = True
    if not skip_case_count_filter:
        criteria &= (ht.n_cases_by_pop >= min_cases)

    if single_sex_only:
        prop_female = ht.n_cases_females / (ht.n_cases_males + ht.n_cases_females)
        criteria &= ((prop_female <= 0.1) | (prop_female >= 0.9))

    ht = ht.filter(criteria).key_by()

    if sex_stratified:
        ht_sex_specific = ht.annotate(pheno_sex='males').union(ht.annotate(pheno_sex='females'))
        if sex_stratified == 'all':
            ht = ht.union(ht_sex_specific)
        else:
            ht = ht_sex_specific

    output = set([tuple(x[field] for field in PHENO_KEY_FIELDS) for x in ht.select(*PHENO_KEY_FIELDS).collect()])
    if specific_phenos:
        specific_phenos = specific_phenos.split(',')
        output = [x for x in output if all(map(lambda y: y is not None, x)) and any([re.match(pcd, '-'.join(x)) for pcd in specific_phenos])]
    if limit:
        output = set(sorted(output)[:limit])

    pheno_key_dict = [dict(zip(PHENO_KEY_FIELDS, x)) for x in output]
    return pheno_key_dict


def custom_summarize_data(suffix, overwrite):
    mt = hl.read_matrix_table(get_custom_ukb_pheno_mt_path(suffix))
    ht = mt.group_rows_by('pop').aggregate(
        stats=hl.agg.stats(mt.both_sexes),
        n_cases_by_pop=hl.cond(hl.set({'continuous', 'biomarkers'}).contains(mt.trait_type),
                               hl.agg.count_where(hl.is_defined(mt.both_sexes)),
                               hl.int64(hl.agg.sum(mt.both_sexes)))
    ).entries()
    ht = ht.key_by('pop', *PHENO_KEY_FIELDS)
    ht = ht.checkpoint(get_custom_phenotype_summary_path(suffix, 'full'), overwrite=overwrite, _read_if_exists=not overwrite)
    ht.flatten().export(get_custom_phenotype_summary_path(suffix, 'full', 'tsv'))


def custom_add_description(mt):
    ht = hl.import_table(DESCRIPTION_PATH).key_by('phenotype_id')
    mt = mt.annotate_cols(description = ht[mt.phesant_pheno].description)
    mt = mt.annotate_cols(description = hl.if_else(hl.is_missing(mt.description), "", mt.description))
    return mt


def custom_load_custom_pheno(data_path, trait_type, modifier, source, sex: str = 'both_sexes', extension: str = 'txt', sample_col='s'):
    print(f'Loading {data_path}...')
    if extension == 'ht':
        ht = hl.read_table(data_path)
    else:
        ht = hl.import_table(data_path, impute=True)
        ht = ht.annotate(**{sample_col: hl.tstr(ht[sample_col])})
        ht = ht.key_by(userId=ht[sample_col])
        if sample_col != 'userId':
            ht = ht.drop(sample_col)
        if trait_type == 'categorical':
            ht = ht.annotate(**{x: hl.bool(ht[x]) for x in list(ht.row_value)})

    mt = pheno_ht_to_mt(ht, trait_type, rekey=False).annotate_cols(data_type=trait_type)
    mt = custom_add_description(mt)
    mt = mt.key_cols_by(trait_type=trait_type, phenocode=mt.phesant_pheno, pheno_sex=sex, coding=NULL_STR_KEY,
                        modifier=modifier).drop('phesant_pheno')
    mt = mt.annotate_cols(category=source)
    return mt


def produce_custom_phenotype_mt(data_path, extn, suffix, trait_type, modifier, source, sample_col='s', append=True, overwrite=False):
    curdate = date.today().strftime("%y%m%d")
    mt = custom_load_custom_pheno(data_path, trait_type=trait_type, modifier=modifier, 
                                  source=source, sample_col=sample_col,
                                  extension=extn
                                  ).checkpoint(get_custom_munged_pheno_path(suffix), args.overwrite)
    cov_ht = get_covariates(hl.int32).persist()
    mt = combine_pheno_files_multi_sex_legacy({'custom': mt}, cov_ht)

    mt.group_rows_by('pop').aggregate(
        n_cases=hl.agg.count_where(mt.both_sexes == 1.0),
        n_controls=hl.agg.count_where(mt.both_sexes == 0.0)
    ).entries().drop(*[x for x in PHENO_COLUMN_FIELDS if x != 'description']).show(100, width=180)
    
    mt_path = get_custom_ukb_pheno_mt_path(suffix)
    if append and hl.hadoop_exists(f'{mt_path}/_SUCCESS'):
        original_mt = hl.read_matrix_table(mt_path)
        original_mt = original_mt.checkpoint(get_custom_ukb_pheno_mt_path(f'{suffix}_before_{curdate}'), overwrite=overwrite)
        original_mt.cols().export(get_custom_phenotype_summary_backup_path(suffix, curdate))
        original_mt.union_cols(mt, row_join_type='outer').write(mt_path, overwrite=overwrite)
    else:
        mt.write(mt_path, overwrite=overwrite)
    custom_summarize_data(suffix, overwrite=overwrite)


def export_pheno_custom(p: Batch, output_path: str, pheno_keys, module: str, mt_loading_function: str,
                        docker_image: str, proportion_single_sex: float = 0.1, n_threads: int = 8, storage: str = '500Mi',
                        additional_args: str = ''):
    extract_task: Job = p.new_job(name='extract_pheno', attributes=copy.deepcopy(pheno_keys))
    extract_task.image(docker_image).cpu(n_threads).storage(storage)
    pheno_dict_opts = ' '.join([f"--{k} {shq(v)}" for k, v in pheno_keys.items()])
    python_command = f"""set -o pipefail; python3 {SCRIPT_DIR}/export_pheno.py
    --load_module {module} --load_mt_function {mt_loading_function}
    {pheno_dict_opts} {"--binary_trait" if saige_pheno_types.get(pheno_keys['trait_type']) != 'quantitative' else ""}
    --proportion_single_sex {proportion_single_sex}
    {"--additional_args " + additional_args if additional_args else ''}
    --output_file {extract_task.out}
    --n_threads {n_threads} | tee {extract_task.stdout}
    ; """.replace('\n', ' ')
    activate_service_account(extract_task)
    extract_task.command(python_command)

    p.write_output(extract_task.out, output_path)
    p.write_output(extract_task.stdout, f'{output_path}.log')
    return extract_task


def main(args):
    hl.init(log='/tmp/saige_temp_hail.log')

    # num_pcs = 20
    num_pcs = 10
    start_time = time.time()
    if args.include_base_covariates:
        basic_covars = ['sex', 'age', 'age2', 'age_sex', 'age2_sex']
    else:
        basic_covars = []
    covariates = ','.join(basic_covars + [f'PC{x}' for x in range(1, num_pcs + 1)])
    n_threads = 8
    analysis_type = "variant"
    chromosomes = list(map(str, range(1, 23))) + ['X']
    reference = 'GRCh37'
    chrom_lengths = hl.get_reference(reference).lengths
    iteration = 1
    pops = args.pops.split(',') if args.pops else POPS

    # if args.local_test:
    #     backend = hb.LocalBackend(gsa_key_file='/Users/konradk/.hail/ukb-diverse-pops.json')
    # else:

    # ensure the data exist
    if not hl.hadoop_exists(f'{get_custom_ukb_pheno_mt_path(args.suffix)}/_SUCCESS') or (args.overwrite_pheno_data or args.append):
        produce_custom_phenotype_mt(args.data_path, args.data_extn, args.suffix, 
                                    trait_type=args.trait_type, modifier=args.modifier,
                                    source=args.source, sample_col=args.sample_col, 
                                    append=args.append, overwrite=args.overwrite_pheno_data)

    backend = hb.ServiceBackend(billing_project='ukb_diverse_pops',
                                bucket=temp_bucket.split('gs://', 1)[-1])
    for pop in pops:
        p = hb.Batch(name=f'saige_pan_ancestry_{pop}', backend=backend, default_image=SAIGE_DOCKER_IMAGE,
                     default_storage='500Mi', default_cpu=n_threads)
        window = '1e7' if pop == 'EUR' else '1e6'
        logger.info(f'Setting up {pop}...')
        chunk_size = int(5e6) if pop != 'EUR' else int(1e6)
        phenos_to_run = custom_get_phenos_to_run(pop, suffix=args.suffix, limit=int(args.local_test),
                                                 specific_phenos=args.phenos,
                                                 skip_case_count_filter=args.skip_case_count_filter,
                                                 sex_stratified=args.sex_stratified)
        logger.info(f'Got {len(phenos_to_run)} phenotypes...')
        if len(phenos_to_run) <= 20:
            logger.info(phenos_to_run)

        pheno_export_dir = f'{pheno_folder}/exported/{args.suffix}/{pop}'
        phenos_already_exported = {}
        if not args.overwrite_pheno_data and hl.hadoop_exists(pheno_export_dir):
            phenos_already_exported = {x['path'] for x in hl.hadoop_ls(pheno_export_dir)}
        pheno_exports = {}

        for pheno_key_dict in phenos_to_run:
            pheno_export_path = get_pheno_output_path(pheno_export_dir, pheno_key_dict, legacy=False)
            if not args.overwrite_pheno_data and pheno_export_path in phenos_already_exported:
                pheno_file = p.read_input(pheno_export_path)
            else:
                # NOTE need to update module and function name
                pheno_task = export_pheno_custom(p, pheno_export_path, pheno_key_dict, 'ukbb_pan_ancestry', 'get_custom_ukb_pheno_mt', 
                                                 PHENO_DOCKER_IMAGE, additional_args=pop, n_threads=n_threads, proportion_single_sex=0)
                pheno_task.attributes.update({'pop': pop})
                pheno_file = pheno_task.out
            pheno_exports[stringify_pheno_key_dict(pheno_key_dict)] = pheno_file
        completed = Counter([isinstance(x, InputResourceFile) for x in pheno_exports.values()])
        logger.info(f'Exporting {completed[False]} phenos (already found {completed[True]})...')

        overwrite_null_models = args.create_null_models
        null_model_dir = f'{root}/null_glmm/{args.suffix}/{pop}'
        null_models_already_created = {}
        if not overwrite_null_models and hl.hadoop_exists(null_model_dir):
            null_models_already_created = {x['path'] for x in hl.hadoop_ls(null_model_dir)}
        null_models = {}

        for pheno_key_dict in phenos_to_run:
            null_glmm_root = get_pheno_output_path(null_model_dir, pheno_key_dict, '', legacy=False)
            model_file_path = f'{null_glmm_root}.rda'
            variance_ratio_file_path = f'{null_glmm_root}.{analysis_type}.varianceRatio.txt'

            if not overwrite_null_models and model_file_path in null_models_already_created and \
                    variance_ratio_file_path in null_models_already_created:
                model_file = p.read_input(model_file_path)
                variance_ratio_file = p.read_input(variance_ratio_file_path)
            else:
                if args.skip_any_null_models: continue
                fit_null_task = fit_null_glmm(p, null_glmm_root, pheno_exports[stringify_pheno_key_dict(pheno_key_dict)],
                                              pheno_key_dict['trait_type'], covariates,
                                              get_ukb_grm_plink_path(pop, iteration, window), SAIGE_DOCKER_IMAGE,
                                              inv_normalize=False, n_threads=n_threads, min_covariate_count=1,
                                              non_pre_emptible=args.non_pre_emptible, storage='100Gi')
                fit_null_task.attributes.update({'pop': pop})
                fit_null_task.attributes.update(copy.deepcopy(pheno_key_dict))
                model_file = fit_null_task.null_glmm.rda
                variance_ratio_file = fit_null_task.null_glmm[f'{analysis_type}.varianceRatio.txt']
            null_models[stringify_pheno_key_dict(pheno_key_dict)] = (model_file, variance_ratio_file)

        completed = Counter([isinstance(x[0], InputResourceFile) for x in null_models.values()])
        logger.info(f'Running {completed[False]} null models (already found {completed[True]})...')

        use_bgen = True
        vcf_dir = f'{root_vcf}/vcf/{pop}'
        test_extension = 'bgen' if use_bgen else 'vcf.gz'
        overwrite_vcfs = args.create_vcfs
        vcfs_already_created = {}
        if not overwrite_vcfs and hl.hadoop_exists(vcf_dir):
            vcfs_already_created = {x['path'] for x in hl.hadoop_ls(vcf_dir)}
        # logger.info(f'Found {len(vcfs_already_created)} VCFs in directory...')
        vcfs = {}
        for chromosome in chromosomes:
            chrom_length = chrom_lengths[chromosome]
            for start_pos in range(1, chrom_length, chunk_size):
                end_pos = chrom_length if start_pos + chunk_size > chrom_length else (start_pos + chunk_size)
                interval = f'{chromosome}:{start_pos}-{end_pos}'
                vcf_root = f'{vcf_dir}/variants_{chromosome}_{str(start_pos).zfill(9)}'
                if f'{vcf_root}.{test_extension}' in vcfs_already_created:
                    if use_bgen:
                        vcf_file = p.read_input_group(**{'bgen': f'{vcf_root}.bgen',
                                                         'bgen.bgi': f'{vcf_root}.bgen.bgi',
                                                         'sample': f'{vcf_root}.sample'})
                    else:
                        vcf_file = p.read_input_group(**{'vcf.gz': f'{vcf_root}.vcf.gz',
                                                         'vcf.gz.tbi': f'{vcf_root}.vcf.gz.tbi'})
                else:
                    vcf_task = extract_vcf_from_mt(p, vcf_root, HAIL_DOCKER_IMAGE, 'ukbb_pan_ancestry', adj=False,
                                                   additional_args=f'{chromosome},{pop}', input_dosage=True,
                                                   reference=reference, interval=interval, export_bgen=use_bgen,
                                                   n_threads=n_threads)
                    vcf_task.attributes['pop'] = pop
                    vcf_file = vcf_task.out
                vcfs[interval] = vcf_file
                if args.local_test:
                    break
            if args.local_test:
                break

        completed = Counter([type(x) == InputResourceFile for x in vcfs.values()])
        logger.info(f'Creating {completed[False]} VCFs (already found {completed[True]})...')

        result_dir = f'{root}/result/{args.suffix}/{pop}'
        overwrite_results = args.overwrite_results
        log_pvalue = True
        for i, pheno_key_dict in enumerate(phenos_to_run):
            if stringify_pheno_key_dict(pheno_key_dict) not in null_models: continue
            model_file, variance_ratio_file = null_models[stringify_pheno_key_dict(pheno_key_dict)]

            if not i % 10:
                n_jobs = dict(Counter(map(lambda x: x.name, p.select_jobs("")))).get("run_saige", 0)
                logger.info(f'Read {i} phenotypes ({n_jobs} new to run so far)...')

            pheno_results_dir = get_pheno_output_path(result_dir, pheno_key_dict, '', legacy=False)
            results_already_created = {}

            if not overwrite_results and not args.skip_saige and hl.hadoop_exists(pheno_results_dir):
                results_already_created = {x['path'] for x in hl.hadoop_ls(pheno_results_dir)}

            saige_tasks = []
            for chromosome in chromosomes:
                if args.skip_saige: break
                chrom_length = chrom_lengths[chromosome]
                for start_pos in range(1, chrom_length, chunk_size):
                    end_pos = chrom_length if start_pos + chunk_size > chrom_length else (start_pos + chunk_size)
                    interval = f'{chromosome}:{start_pos}-{end_pos}'
                    vcf_file = vcfs[interval]
                    results_path = get_results_prefix(pheno_results_dir, pheno_key_dict, chromosome, start_pos, legacy=False)
                    if overwrite_results or f'{results_path}.single_variant.txt' not in results_already_created:
                        samples_file = p.read_input(get_ukb_samples_file_path(pop, iteration))
                        saige_task = run_saige(p, results_path, model_file, variance_ratio_file, vcf_file, samples_file,
                                               SAIGE_DOCKER_IMAGE, trait_type=pheno_key_dict['trait_type'], use_bgen=use_bgen,
                                               chrom=chromosome, log_pvalue=log_pvalue)
                        saige_task.attributes.update({'interval': interval, 'pop': pop})
                        saige_task.attributes.update(copy.deepcopy(pheno_key_dict))
                        saige_tasks.append(saige_task)
                    if args.local_test:
                        break
                if args.local_test:
                    break

            res_tasks = []
            if overwrite_results or args.overwrite_hail_results or \
                    f'{pheno_results_dir}/variant_results.ht' not in results_already_created or \
                    not hl.hadoop_exists(f'{pheno_results_dir}/variant_results.ht/_SUCCESS'):
                null_glmm_root = get_pheno_output_path(null_model_dir, pheno_key_dict, f'.{analysis_type}.log',
                                                       legacy=False)

                prefix = get_results_prefix(pheno_results_dir, pheno_key_dict,
                                            f'{"chr" if reference == "GRCh38" else ""}{{chrom}}', 1,
                                            legacy=False)
                saige_log = f'{prefix}.{analysis_type}.log'

                load_task = load_results_into_hail(p, pheno_results_dir, pheno_key_dict,
                                                   saige_tasks, get_ukb_vep_path(), HAIL_DOCKER_IMAGE,
                                                   saige_log=saige_log, analysis_type=analysis_type,
                                                   n_threads=n_threads, null_glmm_log=null_glmm_root,
                                                   reference=reference, legacy_annotations=True,
                                                   log_pvalue=log_pvalue)
                load_task.attributes['pop'] = pop
                res_tasks.append(load_task)
                qq_export, qq_plot = qq_plot_results(p, pheno_results_dir, res_tasks, HAIL_DOCKER_IMAGE, QQ_DOCKER_IMAGE, n_threads=n_threads)
                qq_export.attributes.update({'pop': pop})
                qq_export.attributes.update(copy.deepcopy(pheno_key_dict))
                qq_plot.attributes.update({'pop': pop})
                qq_plot.attributes.update(copy.deepcopy(pheno_key_dict))

        logger.info(f'Setup took: {time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))}')
        logger.info(f'Submitting: {get_tasks_from_pipeline(p)}')
        logger.info(f"Total size: {sum([len(x._pretty()) for x in p.select_jobs('')])}")
        p.run(dry_run=args.dry_run, wait=False, delete_scratch_on_exit=False)
        logger.info(f'Finished: {get_tasks_from_pipeline(p)}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--single_sex_only', help='Run only single sex phenotypes (experimental)', action='store_true')
    parser.add_argument('--sex_stratified', help='Run these phenotypes in a sex-stratified fashion (experimental)', choices=(None, 'all', 'only'))
    parser.add_argument('--skip_any_null_models', help='Skip running SAIGE null models', action='store_true')
    parser.add_argument('--skip_saige', help='Skip running SAIGE tests', action='store_true')
    parser.add_argument('--create_null_models', help='Force creation of null models', action='store_true')
    parser.add_argument('--create_vcfs', help='Force creation of VCFs', action='store_true')
    parser.add_argument('--overwrite_pheno_data', help='Overwrite phenotype munged data and exports', action='store_true')
    parser.add_argument('--overwrite_results', help='Force run of SAIGE tests', action='store_true')
    parser.add_argument('--overwrite_hail_results', help='Force run of results loading', action='store_true')
    parser.add_argument('--local_test', help='Local test of pipeline', action='store_true')
    parser.add_argument('--non_pre_emptible', help='Local test of pipeline', action='store_true')
    parser.add_argument('--skip_case_count_filter', help='Skip running SAIGE tests', action='store_true')
    parser.add_argument('--phenos', help='Comma-separated list of trait_type-phenocode-pheno_sex-coding-modifier regexes '
                                         '(e.g. continuous-50-both_sexes--,icd10-E1.*,brain_mri-.* )')
    parser.add_argument('--pops', help='comma-searated list')
    parser.add_argument('--dry_run', help='Dry run only', action='store_true')
    parser.add_argument('--include_base_covariates', help='If true, will include the usual age, sex covariates.', action='store_true')
    parser.add_argument('--data_path', required=True, type=str, help='Path to phenotype data.')
    parser.add_argument('--sample_col', default='s', type=str, help='Sample identifier column.')
    parser.add_argument('--data_extn', required=True, type=str, help='Type of the phenotype data (ht or tsv).')
    parser.add_argument('--suffix', required=True, type=str, help='Analysis suffix.')
    parser.add_argument('--trait_type', required=True, type=str, help='Trait type; also used for munging.')
    parser.add_argument('--modifier', required=True, type=str, help='Global modifier for this dataset used for phenotype munging.')
    parser.add_argument('--source', required=True, type=str, help='Data source for tagging, used for pheontype munging.')
    parser.add_argument('--append', action='store_true', help='If enabled, will attempt to augment the pheno MT with new phenotypes.')
    args = parser.parse_args()

    main(args)
