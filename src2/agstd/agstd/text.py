import logging
log = logging.getLogger(__name__)

import numpy as np

class TXTParser(object):
    """
    This Object facilitate the parsing of a standard table from a text file
    loading the data into a numpy array based on the given description and
    numpy data type.

    Optionally we can provide our own object and line delimiter (nearly CSV)

    inputs :
        description - A list of 4-tuple representing :
            1 - column name
            2 - a callable to parse the column argument
            3 - A callable to create a beautiful output from the python structure
            4 - The name of the column associated in the numpy dtype

    """
    def __init__(self, description, dtype, object_delimiter = " ", line_delimiter = "\n"):
        self.description = description
        self.line_delimiter = line_delimiter
        self.object_delimiter = object_delimiter
        self.ids = {}
        self.dtype = dtype
        for i, (name, intype, outype, ident, cname) in enumerate(description):
            self.ids[name] = i

    def loads(self, txt):
        """
        This method load a table from a string and return the numpy array
        created.
        """
        lines = txt.split(self.line_delimiter)
        header = lines[0]
        olist = []
        for l in lines[1:]:
            if l == '':
                continue
            desc = l.split()
            o = [intype(desc[i]) for i, (name, intype, outype, ident, csize) in enumerate(self.description)]
            olist.append(o)

        # Creating the actual numpy ary
        ary = np.empty(len(olist), dtype = self.dtype)

        for i in range(len(olist)):
            ary[i] = tuple(olist[i])

        return ary

    def dumps(self, olst, include_header = True):
        """
        This method dumps the content of a table to a human readeable txt
        formated  string.
        """
        if include_header:
            header = []
            for name, intype, outype, ident, csize in self.description:
                header.extend([self.object_delimiter * (csize - len(name)), name])
        else:
            header = ''

        for obj in olst:
            line = [self.line_delimiter]
            for (name, intype, outype, ident, csize) in self.description:
                value = outype(obj[ident])
                line.extend([self.object_delimiter * (csize - len(value)), value])
            header.extend(line)

        return ''.join(header)

    def dump(self, f, *args, **kw):
        f = open(f, 'w') if isinstance(f, str) else f
        f.write(self.dumps(*args, **kw))


    def load(self, f, *args, **kw):
        f = open(f) if isinstance(f, str) else f
        return self.loads(f, *args, **kw)
