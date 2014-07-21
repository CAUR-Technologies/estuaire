def dictify(el):
    """
    This method take a tree element and converts it into a python dictionary.

    WARNING : This method take into account that there is no property with
    the same name as the children nodes. A dictionary cannot represent the
    full extent of the XML Standard ... think about it.
    """
    d = {};
    d.update(el.attrib)
    for k, l in el.__dict__.iteritems():
        d[k] = [dictify(e) for e in l]
    return d
