class InputInChI:

    def input_inchi(get_ids_function, config):
        print("\nEnter first InChI:")
        inchi1 = input().strip()

        print("\nEnter second InChI:")
        inchi2 = input().strip()

        result = get_ids_function(inchi1, inchi2, config)

        print("RESULT:")
        for layer, value in result.items():
            print(f"{layer.name}: {value}")

        return result