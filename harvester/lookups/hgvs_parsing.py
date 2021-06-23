import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.normalizer
import hgvs.parser

hdp = hgvs.dataproviders.uta.connect()
hgnorm = hgvs.normalizer.Normalizer(hdp)
hgvsparser = hgvs.parser.Parser()
am = {
    'GRCh37': hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh37', normalize=True),
    'GRCh38': hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh38', normalize=True)
}
