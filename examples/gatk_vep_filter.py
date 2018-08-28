vep_csq_desc_prefix = "Consequence annotations from Ensembl VEP. Format: "
csq_format_keys = None

important_genes = set(['MSH2', 'MSH2', 'MSH6', 'MLH1', 'PMS2', 'MSH4', 'MUTYH', 'EXO1', 'PMS1', 'MSH3', 'MSH5', 'PMS2P3', 'POLE', 'MLH3', 'FAN1', 'POLD1', 'NTHL1', 'RNF43', 'APC', 'BRCA1', 'BRCA2', 'AXIN2'])


def parse_gt(text):
    fields = text.split(':')
    return fields[0], fields[1]

def filter(record, metadata, sample_ids):
    global vep_csq_desc_prefix
    global csq_format_keys
    global important_genes
    global parse_gt
    genotypes = [parse_gt(gt) for gt in record.genotypes]
    num_carriers = sum([1 for gt in genotypes if '1' in gt[0]]) 
    genotypes_strs = ["{}:{}".format(gt[0], gt[1]) for gt in genotypes]
    if num_carriers <= 6:
        if csq_format_keys is None:
            csq_format_keys = (metadata['INFO']['CSQ']['Description'][len(vep_csq_desc_prefix):]).split('|')
        picked_csq = {} 
        if 'CSQ' in record.info:
            csqs = record.info['CSQ']
            if type(csqs) != list:
                csqs = [csqs]
            for csq in csqs:
                csq_values = csq.split('|')
                csq_dict = dict(zip(csq_format_keys, csq_values))
                if csq_dict['PICK'] == '1':
                    picked_csq = csq_dict
                    break
        consequence = picked_csq.get('Consequence', '.')
        impact = picked_csq.get('IMPACT', '.')
        gene_symbol = picked_csq.get('SYMBOL', '.')
        exon = picked_csq.get('EXON', '.')
        polyphen = picked_csq.get('PolyPhen', '.')
        sift = picked_csq.get('SIFT', '.')
        c_change = picked_csq.get('HGVSc', '.')
        p_change = picked_csq.get('HGVSp', '.')
        gnomad_af = picked_csq.get('gnomAD_AF', '.')
        gnomad_af_float = 0.0
        important_gene = '1' if gene_symbol in important_genes else '0'
        try:
            gnomad_af_float = float(gnomad_af)
        except:
            pass
        if gnomad_af_float <= 0.05 and impact != 'LOW' and consequence not in [".", 'intergenic_variant', 'intron_variant', 'downstream_gene_variant', 'upstream_gene_variant', 'intron_variant&non_coding_transcript_variant', 'non_coding_transcript_exon_variant']:
            row = [record.chrom, record.pos, record.ref, record.alt, record.qual, consequence, impact, gene_symbol, important_gene, exon, c_change, p_change, polyphen, sift, gnomad_af, str(num_carriers)] + genotypes_strs
            yield '\t'.join(row)


def begin(_metadata, sample_ids):
    header = ["chrom", "pos", "ref", "alt", "qual", "consequence", "impact", "gene", "CRC gene", "exon", "HGVSc", "HGVSp", "polyphen", "sift", "gnomad_af", "num_carriers"] + sample_ids
    yield "\t".join(header)


def end(_metadata, _sample_ids):
    return iter(())
