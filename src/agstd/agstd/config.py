import ConfigParser

description = {'good' : {'this' : int}}

class TypedConfigParser(ConfigParser.ConfigParser):
    """
    This class is simply intended to make explicit the type of the
    expected configuration item.
    """
    description = {}
    def __init__(self, description, *args, **kw):
        ConfigParser.ConfigParser.__init__(self, *args, **kw)
        self.description = description

    def get(self, section, item, *args, **kw):
        value = ConfigParser.ConfigParser.get(self, section, item, *args, **kw)
        if section in self.description:
            section_description = self.description[section]
            if item in section_description:
                value = section_description[item](value)
        return value


