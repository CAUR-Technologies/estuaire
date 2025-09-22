import os
import pickle
import weakref


class DirectoryDB(object):
    @staticmethod
    def join_keys(args):
        if hasattr(args, "__iter__"):
            key = ""
            for item in args:
                key = f"{key}/{item}"
            key = key[1:]
        else:
            key = args
        return key

    def __init__(self, root, protocol=pickle.HIGHEST_PROTOCOL):
        self.root = os.path.expanduser(root)
        if not os.path.exists(self.root):
            os.makedirs(self.root)
        self.protocol = protocol
        self.cache = {}

    def __len__(self):
        return len(self.keys())

    def __iter__(self):
        return iter(self.keys())

    def __getitem__(self, keys):
        key = self.join_keys(keys)
        if key in self.cache:
            item = self.cache[key]()
            if item is not None:
                return item

        keypath = os.path.join(self.root, key)
        if not os.path.exists(keypath):
            raise KeyError(key)

        if os.path.isdir(keypath):
            return DirectoryDB(keypath)

        with open(keypath, "rb") as fdesc:
            value = pickle.load(fdesc)

        try:
            self.cache[key] = weakref.ref(value)
        except TypeError:
            pass
        return value

    def keys(self):
        return os.listdir(self.root)

    def __setitem__(self, keys, value):
        self.set(keys, value)

    def set(self, keys, value, protocol=None):
        key = self.join_keys(keys)
        if protocol is None:
            protocol = self.protocol
        split_path = os.path.split(key)
        dirpath = os.path.join(self.root, *split_path[:-1])
        keypath = os.path.join(self.root, *split_path)
        if not os.path.exists(dirpath):
            os.makedirs(dirpath)
        with open(keypath, "wb") as fdesc:
            pickle.dump(value, fdesc, protocol=protocol)
        try:
            self.cache[key] = weakref.ref(value)
        except TypeError:
            pass

    def __delitem__(self, key):
        if key in self.cache:
            del self.cache[key]
        keypath = os.path.join(self.root, key)
        if not os.path.exists(keypath):
            raise KeyError(key)
        os.remove(keypath)
