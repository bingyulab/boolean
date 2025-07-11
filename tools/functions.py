
    # Convert R list to Python dict
    def rlist_to_pydict(rlist):
        py_dict = {}
        for name in rlist.names:
            value = rlist.rx2(name)
            # Convert R vectors to Python lists
            if hasattr(value, 'tolist'):
                py_dict[name] = value.tolist()
            else:
                py_dict[name] = value
        return py_dict