import hail as hl

bucket = 'gs://ukb-diverse-pops'
REFERENCE_GENOME = 'GRCh37'
CHROMOSOMES = list(map(str, range(1, 23))) + ['X', 'XY']
POPS = ['AFR', 'EAS', 'CSA', 'MID', 'EUR', 'AMR']


def get_hq_samples():
    ht = hl.import_table(f'{bucket}/misc/ukb31063_samples_qc_FULL.txt', no_header=True)
    return ht.key_by(s=ht.f0).drop('f0')


def get_covariates_tsv_path():
    return f'{bucket}/pca/non_eur_within_pop_pc_covs.txt.bgz'


def get_covariates(key_type = hl.tstr):
    return hl.import_table(get_covariates_tsv_path(), key='s', impute=True, types={'s': key_type})


def get_ukb_meta_pop_tsv_path():
    return f'{bucket}/pca/globalref_ukbb_pca_pops_rf_50.txt.bgz'


def get_ukb_meta(key_type = hl.tstr):
    return hl.import_table(get_ukb_meta_pop_tsv_path(), key='s', impute=True, types={'s': key_type})


def get_ukb_pheno_mt(pop: str = 'all'):
    from .phenotypes import get_ukb_pheno_mt_path
    mt = hl.read_matrix_table(get_ukb_pheno_mt_path())
    mt = mt.annotate_rows(**get_ukb_meta(key_type=hl.tint32)[mt.row_key])
    if pop != 'all':
        mt = mt.filter_rows(mt.pop == pop)
    return mt

