
class Harvester(object):
    def __init__(self):
        pass

    def harvest(self, genes):
        raise NotImplementedError()

    def convert(self, gene_data):
        raise NotImplementedError()

    def harvest_and_convert(self, genes):
        raise NotImplementedError()

    def test(self):
        print "No implemented test for harvester %s" % str(self)
        return

    def __str__(self):
        return "Harvester"
