import hail as hl
import argparse
import os
import sys
import tempfile
import warnings
import importlib


def gt_to_gp(mt, location: str = 'GP'):
    return mt.annotate_entries(**{location: hl.or_missing(
        hl.is_defined(mt.GT),
        hl.map(lambda i: hl.cond(mt.GT.unphased_diploid_gt_index() == i, 1.0, 0.0),
               hl.range(0, hl.triangle(hl.len(mt.alleles)))))})


def impute_missing_gp(mt, location: str = 'GP', mean_impute: bool = True):
    mt = mt.annotate_entries(_gp = mt[location])
    if mean_impute:
        mt = mt.annotate_rows(_mean_gp=hl.agg.array_agg(lambda x: hl.agg.mean(x), mt._gp))
        gp_expr = mt._mean_gp
    else:
        gp_expr = [1.0, 0.0, 0.0]
    return mt.annotate_entries(**{location: hl.or_else(mt._gp, gp_expr)}).drop('_gp')


def main(args):
    hl.init(master=f'local[{args.n_threads}]',
            log=hl.utils.timestamp_path(os.path.join(tempfile.gettempdir(), 'extract_vcf'), suffix='.log'),
            default_reference=args.reference)

    sys.path.append('/')
    add_args = []
    if args.additional_args is not None:
        add_args = args.additional_args.split(',')
    load_module = importlib.import_module(args.load_module)
    mt = getattr(load_module, args.load_mt_function)(*add_args, **{'entry_fields': ('GP','GT',)})

    if args.gene_map_ht_path is None:
        interval = [hl.parse_locus_interval(args.interval)]
    else:
        gene_ht = hl.read_table(args.gene_map_ht_path)
        if args.gene is not None:
            gene_ht = gene_ht.filter(gene_ht.gene_symbol == args.gene)
            interval = gene_ht.aggregate(hl.agg.take(gene_ht.interval, 1), _localize=False)
        else:
            gene_ht = hl.filter_intervals(gene_ht, [hl.parse_locus_interval(args.gene_ht_interval)])
            interval = [hl.parse_locus_interval(args.interval)]

        gene_ht = gene_ht.filter(hl.set(args.groups.split(',')).contains(gene_ht.annotation))
        if args.common_variants_only:
            gene_ht = gene_ht.filter(gene_ht.common_variant)
        gene_ht.select(group=gene_ht.gene_id + '_' + gene_ht.gene_symbol + '_' + gene_ht.annotation +
                       hl.if_else(gene_ht.common_variant, '_' + gene_ht.variants[0], ''),
                       variant=hl.delimit(gene_ht.variants, '\t')
                       ).key_by().drop('start').export(args.group_output_file, header=False)
        # TODO: possible minor optimization: filter output VCF to only variants in `gene_ht.variants`

    if args.group_file_only:
        return

    if not args.no_adj:
        mt = mt.filter_entries(mt.adj)

    mt = hl.filter_intervals(mt, interval)

    if not args.input_bgen or args.encoding == 'recessive':
        mt = mt.select_entries('GT')
    else:
        mt = mt.select_entries('GP')
    if not args.input_bgen:
        mt = mt.filter_rows(hl.agg.count_where(mt.GT.is_non_ref()) > 0)
    mt = mt.annotate_rows(rsid=mt.locus.contig + ':' + hl.str(mt.locus.position) + '_' + mt.alleles[0] + '/' + mt.alleles[1])

    if args.callrate_filter:
        mt = mt.filter_rows(hl.agg.fraction(hl.is_defined(mt.GT)) >= args.callrate_filter)

    if args.export_bgen:
        if not args.input_bgen or args.encoding == 'recessive':
            mt = mt.annotate_entries(GT=hl.if_else(mt.GT.is_haploid(), hl.call(mt.GT[0], mt.GT[0]), mt.GT))
            if args.encoding == 'recessive':
                mt = mt.annotate_entries(GT = hl.if_else(mt.GT.is_het(), hl.call(0, 0, phased=False), mt.GT))
            mt = gt_to_gp(mt)
            mt = impute_missing_gp(mt, mean_impute=args.mean_impute_missing)
        
        if mt.count_rows() == 0:
            with open(args.output_file + '.bgen', 'w') as f:
                f.write('Failed due to no variants.')
            with open(args.output_file + '.sample', 'w') as f:
                f.write('Failed due to no variants.')
            with open(args.output_file + '.bgen.bgi', 'w') as f:
                f.write('Failed due to no variants.')
            with open(args.output_file + '.failure', 'w') as f:
                f.write('Failed due to no variants.')
        else:
            hl.export_bgen(mt, args.output_file, gp=mt.GP, varid=mt.rsid)
    else:
        mt = mt.annotate_entries(GT=hl.or_else(mt.GT, hl.call(0, 0)))
        # Note: no mean-imputation for VCF
        hl.export_vcf(mt, args.output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--load_module', help='Module to load with helper functions', required=True)
    parser.add_argument('--load_mt_function', help='Function in module to load analysis-ready MatrixTable', default='get_filtered_mt')
    parser.add_argument('--additional_args', help='Comma separated list of arguments to pass to mt_function')

    parser.add_argument('--gene_map_ht_path', help='Path to gene map HT')
    parser.add_argument('--groups', help='Which variant groupings to use')
    parser.add_argument('--gene', help='Gene to export')
    parser.add_argument('--interval', help='Interval to export')
    parser.add_argument('--gene_ht_interval', help='Interval to export')

    parser.add_argument('--input_bgen', help='Input is BGEN instead of MT', action='store_true')
    parser.add_argument('--no_adj', help='Use all genotypes instead of only high-quality ones', action='store_true')
    parser.add_argument('--mean_impute_missing', help='Whether to mean impute missing genotypes (BGEN only) '
                                                      '(default: set to hom ref)', action='store_true')
    parser.add_argument('--export_bgen', help='Export BGEN instead of VCF', action='store_true')
    parser.add_argument('--callrate_filter', help='Impose filter of specified callrate (default: none)', default=0.0, type=float)
    parser.add_argument('--reference', help='Reference genome to use', default='GRCh38', choices=('GRCh37', 'GRCh38'))

    parser.add_argument('--group_output_file', help='Output file for variant groupings')
    parser.add_argument('--group_file_only', help='Only output variant groupings (no VCF or BGEN)', action='store_true')
    parser.add_argument('--common_variants_only', help='Only output common variants', action='store_true')
    parser.add_argument('--output_file', help='Output file')
    parser.add_argument('--n_threads', help='Number of threads', type=int, default=8)
    parser.add_argument('--encoding', help='Additive or recessive encoding', type=str, default='additive')
    args = parser.parse_args()

    if int(args.gene is not None) + int(args.interval is not None) != 1:
        sys.exit('Error: One and only one of --gene or --interval must be specified')

    if args.gene_map_ht_path is None:
        if args.interval is None:
            sys.exit('Error: If --gene_map_ht_path is not specified, --interval must be specified')
    else:
        if args.groups is None or args.group_output_file is None:
            sys.exit('Error: If --gene_map_ht_path is specified, --groups and --group_output_file must be specified')

    main(args)


