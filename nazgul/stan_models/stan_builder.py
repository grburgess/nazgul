import collections

_allowed_block_names = [
    "functions",
    "data",
    "transformed data",
    "parameters",
    "transformed parameters",
    "model",
    "generated quantities",
]


class StanBlock(object):
    def __init__(self, block_name):

        assert block_name in _allowed_block_names, "%s must be one of %s" "".join(
            _allowed_block_names
        )
        self._code_string = f"{block_name} {\n"

        self._generated = False

    def add_include(self, include_file):
        """
        add include files into the block.
        """
        self._code_string += f"#include {include_file}\n"

    def insert_code(self, code):

        self._code_string += f"\n{code}\n"

    def clean(self):

        self._code_string = self._code_string[:-4]
        self._generated = False

    def generate(self):

        self._code_string += "\n}\n\n"
        self._generated = True

    @property
    def code(self):
        if self._generated:

            return self._code_string

        else:

            raise RuntimeError("Block is not generated")


class FunctionsBlock(StanBlock):
    def __init__(self):

        super(FunctionsBlock, self).__init__(block_name="functions")


class DataBlock(StanBlock):
    def __init__(self, data_size="N"):

        super(DataBlock, self).__init__(block_name="data")

        self._data_size = data_size

        # add on the data size
        self._code_string += f" int {data_size};\n" 

    def add_vector_data(self, name, size=None,n_vectors=None stan_type="vector"):

        if size is None:

            size = self._data_size

        if stan_type == "vector":

            if n_vectors is None:
            
                self._code_string += f"  {stan_type}[{size}] {name};\n"

            else:

                self._code_string += f"  {stan_type}[{size}] {name}[n_vectors];\n"

        elif stan_type == "int":

            if n_vectors is None:
            
                self._code_string += f"  {stan_type} {name}[{size}];\n"

            else:
                 self._code_string += f"  {stan_type} {name}[{n_vectors},{size}];\n"
                

                

                
    def add_data(self, name, stan_type="real"):

        self._code_string += f"  {stan_type} {name};\n"


class TransDataBlock(DataBlock):
    def __init__(self):

        super(DataBlock, self).__init__(block_name="transformed data")


class ParametersBlock(StanBlock):
    def __init__(self, data_size="N"):

        super(ParametersBlock, self).__init__(block_name="parameters")

        self._data_size = data_size

    @staticmethod
    def bounds_generator(upper_bound, lower_bound):

        bounds = ""

        if (lower_bound is not None) or (upper_bound is not None):
            bounds += "<"

            if lower_bound is not None:
                bounds += f"lower={lower_bound},"

            if upper_bound is not None:

                bounds += f"upper={upper_bound}"

            else:
                bounds = bounds.replace(",", "")

            bounds += ">"

        return bounds

    def add_vector_parameters(
        self, name, size=None, lower_bound=None, upper_bound=None, stan_type="vector"
    ):

        if size is None:

            size = self._data_size

        bounds = ParametersBlock.bounds_generator(upper_bound, lower_bound)

        self._code_string += f"  {stan_type}{bounds}[{size}] {name};\n"

    def add_parameters(
        self, name, lower_bound=None, upper_bound=None, stan_type="real"
    ):

        bounds = ParametersBlock.bounds_generator(upper_bound, lower_bound)

        self._code_string += f"  {stan_type}{bounds} {name};\n"


class TransParametersBlock(ParametersBlock):
    def __init__(self):

        super(ParametersBlock, self).__init__(block_name="transformed parameters")


class ModelBlock(StanBlock):
    def __init__(self):

        super(ModelBlock, self).__init__(block_name="model")


class GQBlock(StanBlock):
    def __init__(self):

        super(GQBlock, self).__init__(block_name="generated quantities")


class StanGenerator(object):
    def __init__(self, model_name="model", data_size="N"):

        self._file_name = f"{model_name}.stan"
        self._model_name = model_name

        self._blocks = collections.OrderedDict(
            functions=FunctionsBlock(),
            data=DataBlock(data_size=data_size),
            transformed_data=TransDataBlock(),
            parameters=ParametersBlock(data_size=data_size),
            transformed_parameters=TransParametersBlock(),
            model=ModelBlock(),
            generated_quantities=GQBlock(),
        )

    def add_standard_vector_data(self, *data_names):
        """
        add vector data that is the size of the standard data size
        """

        for name in data_names:

            self._blocks["data"].add_vector_data(name)

    def add_vector_data(self, size="M", stan_type="vector", *data_names):
        """
        add vector data that is of size size
        """

        for name in data_names:

            self._blocks["data"].add_vector_data(name, size=size, stan_type=stan_type)

    def add_data(self, *data_names, stan_type="real"):

        for name in data_names:

            self._blocks["data"].add_data(name, stan_type=stan_type)

    def add_standard_vector_parameters(
        self, *parameters_names, lower_bound=None, upper_bound=None
    ):
        """
        add vector parameters that is the size of the standard parameters size
        """

        for name in parameters_names:

            self._blocks["parameters"].add_vector_parameters(
                name, lower_bound=lower_bound, upper_bound=upper_bound
            )

    def add_vector_parameters(
        self, *parameters_names, size="M", lower_bound=None, upper_bound=None
    ):
        """
        add vector parameters that is of size size
        """

        for name in parameters_names:

            self._blocks["parameters"].add_vector_parameters(
                name, size=size, lower_bound=lower_bound, upper_bound=upper_bound
            )

    def add_parameters(
        self, *parameters_names, stan_type="real", lower_bound=None, upper_bound=None
    ):

        for name in parameters_names:

            self._blocks["parameters"].add_parameters(
                name,
                stan_type=stan_type,
                lower_bound=lower_bound,
                upper_bound=upper_bound,
            )

    def write_stan_code(self):

        output_code = ""
        for k, v in self._blocks.items():

            v.generate()
            output_code += v.code
            output_code += "\n"
            v.clean()
        with open(self._file_name, "w") as f:

            f.write(output_code)

    @property
    def blocks(self):

        return self._blocks

    @property
    def data(self):

        self._blocks["data"].generate()
        print(self._blocks["data"].code)
        self._blocks["data"].clean()

    @property
    def parameters(self):

        self._blocks["parameters"].generate()
        print(self._blocks["parameters"].code)
        self._blocks["parameters"].clean()

    @property
    def functions(self):

        self._blocks["functions"].generate()
        print(self._blocks["functions"].code)
        self._blocks["functions"].clean()

    @property
    def generated_quantities(self):

        self._blocks["generated_quantities"].generate()
        print(self._blocks["generated_quantities"].code)
        self._blocks["generated_quantities"].clean()

    @property
    def model(self):

        self._blocks["model"].generate()
        print(self._blocks["model"].code)
        self._blocks["model"].clean()

    @property
    def transformed_data(self):

        self._blocks["transformed_data"].generate()
        print(self._blocks["transformed_data"].code)
        self._blocks["transformed_data"].clean()

    @property
    def transformed_parameters(self):

        self._blocks["transformed_parameters"].generate()
        print(self._blocks["transformed_parameters"].code)
        self._blocks["transformed_parameters"].clean()
