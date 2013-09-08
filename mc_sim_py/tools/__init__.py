class DummyMonomer(object):
    position = (0,0,0)
    type = ''
    neighbours_list = set()
    nb_list = []
    nb_types = []
    nb_types_count = {}

    on_interface = False

    nb_count = 0

    def __init__(self, pos, type, nb):
        self.position = pos
        self.type = type
        self.neighbours_list = set(nb)
        self.nb_types = []
        self.nb_types_count = {}
        self.nb_list = []

    def __str__(self):
        return "%s %d %d %d" % (self.type, self.position[0], self.position[1], self.position[2])

    def __repr__(self):
        return str(self.__str__())
