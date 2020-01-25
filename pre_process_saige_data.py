#!/usr/bin/env python3

__author__ = 'konradk'

from gnomad_hail import *
from ukb_common import *
from ukbb_pan_ancestry import *


def main(args):
    hl.init(default_reference='GRCh37')
    pops = POPS
    pops.remove('EUR')

    if args.create_plink_file:
        for pop in pops:
            call_stats_ht = hl.read_table(get_ukb_af_ht_path(with_x = False))
            mt = get_filtered_mt(pop=pop, imputed=False)
            n_samples = mt.count_cols()
            mt = filter_to_autosomes(mt)
            callstats = call_stats_ht[mt.row_key]
            mt = mt.filter_rows((callstats.an[pop] > 0.95 * n_samples) & (callstats.af[pop] > 0.01))

            mt = mt.checkpoint(get_ukb_grm_mt_path(pop), _read_if_exists=not args.overwrite, overwrite=args.overwrite)
            mt = mt.unfilter_entries()
            ht = hl.ld_prune(mt.GT, r2=0.1)
            ht = ht.checkpoint(get_ukb_grm_pruned_ht_path(pop), _read_if_exists=not args.overwrite, overwrite=args.overwrite)
            mt = mt.filter_rows(hl.is_defined(ht[mt.row_key]))

            if args.overwrite or not hl.hadoop_exists(f'{get_ukb_grm_plink_path(pop)}.bed'):
                hl.export_plink(mt, get_ukb_grm_plink_path(pop))
            print(pop)
            mt = get_filtered_mt(chrom='22', pop=pop)
            if args.overwrite or not hl.hadoop_exists(get_ukb_samples_file_path(pop)):
                with hl.hadoop_open(get_ukb_samples_file_path(pop), 'w') as f:
                    f.write('\n'.join(mt.s.collect()) + '\n')

    if args.vep:
        call_stats_ht = hl.read_table(get_ukb_af_ht_path())
        ht = vep_or_lookup_vep(call_stats_ht)
        ht.write(get_ukb_vep_path(), args.overwrite)

    if args.prepare_genotype_data:
        load_all_mfi_data().write(ukb_imputed_info_ht_path, args.overwrite)

    if args.genotype_summary:
        variants = hl.read_table(ukb_imputed_info_ht_path)
        print(variants.count())
        variants = variants.filter(variants.info > 0.8)
        print(variants.count())
        meta_ht = hl.import_table(get_ukb_meta_pop_tsv_path(), impute=True, types={'s': hl.tstr}, key='s')

        mt = get_ukb_imputed_data('all', variant_list=variants, entry_fields=('dosage', ))
        mt = mt.annotate_cols(**meta_ht[mt.col_key])
        ht = mt.annotate_rows(af=hl.agg.group_by(mt.pop, hl.agg.mean(mt.dosage)),
                              an=hl.agg.group_by(mt.pop, hl.agg.count_where(hl.is_defined(mt.dosage)))).rows()
        ht = ht.checkpoint(get_ukb_af_ht_path(False), args.overwrite, _read_if_exists=not args.overwrite)

        mt = get_ukb_imputed_data('X', variant_list=variants, entry_fields=('dosage', ))
        mt = mt.annotate_cols(**meta_ht[mt.col_key])
        ht_x = mt.annotate_rows(af=hl.agg.group_by(mt.pop, hl.agg.mean(mt.dosage)),
                              an=hl.agg.group_by(mt.pop, hl.agg.count_where(hl.is_defined(mt.dosage)))).rows()
        ht = ht.union(ht_x)
        ht = ht.checkpoint(get_ukb_af_ht_path(), args.overwrite, _read_if_exists=not args.overwrite)

        print(ht.aggregate(hl.struct(
            # hist=hl.agg.hist(hl.sum(ht.an.values()), 0, total_samples, 10),  # No missing data
            # fraction_missingness=hl.agg.fraction(hl.sum(ht.an.values()) < total_samples),
            # number_sites_above_001=hl.agg.count_where(hl.any(lambda x: x / 2 > 0.001, ht.af.values())),
            # number_sites_above_005=hl.agg.count_where(hl.any(lambda x: x / 2 > 0.005, ht.af.values())),
            # number_sites_above_01=hl.agg.count_where(hl.any(lambda x: x / 2 > 0.01, ht.af.values())),
            # number_sites_above_05=hl.agg.count_where(hl.any(lambda x: x / 2 > 0.05, ht.af.values())),
            # number_sites_above_10=hl.agg.count_where(hl.any(lambda x: x / 2 > 0.1, ht.af.values()))
            number_sites_above_mac_20=hl.agg.count_where(hl.any(lambda x: ht.af[x] * ht.an[x] >= 20, hl.literal(POPS))),
            number_run_sites_above_mac_20=hl.agg.sum(hl.sum(hl.map(lambda x: hl.int(ht.af[x] * ht.an[x] >= 20), hl.literal(POPS))))
        )))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--pheno_summary', action='store_true')
    parser.add_argument('--prepare_genotype_data', action='store_true')
    parser.add_argument('--genotype_summary', action='store_true')
    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--create_plink_file', help='Overwrite everything', action='store_true')
    parser.add_argument('--vep', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)

