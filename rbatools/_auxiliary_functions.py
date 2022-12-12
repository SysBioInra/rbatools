# python 2/3 compatibility
from __future__ import division, print_function

import numpy
import difflib
import pandas
import collections
import scipy
import warnings


def return_compartments(model_structure) -> list:
    return(list(model_structure.CompartmentInfo.Elements.keys()))


def return_enzymes(model_structure) -> list:
    return(list(model_structure.EnzymeInfo.Elements.keys()))


def return_proteins(model_structure) -> list:
    return(list(model_structure.ProteinInfo.Elements.keys()))


def return_reactions(model_structure) -> list:
    return(list(model_structure.ReactionInfo.Elements.keys()))


def return_metabolites(model_structure) -> list:
    return(list(model_structure.MetaboliteInfo.Elements.keys()))


def return_processes(model_structure) -> list:
    return(list(model_structure.ProcessInfo.Elements.keys()))


def return_macro_molecules(model_structure) -> list:
    return(list(model_structure.MacromoleculeInfo.Elements.keys()))


def return_targets(model_structure) -> list:
    return(list(model_structure.TargetInfo.Elements.keys()))


def return_modules(model_structure) -> list:
    return(list(model_structure.ModuleInfo.Elements.keys()))


def return_density_constraints(model_structure) -> list:
    return(list(model_structure.DensityConstraintsInfo.Elements.keys()))


def return_process_constraints(model_structure) -> list:
    return(list(model_structure.ProcessConstraintsInfo.Elements.keys()))


def return_enzyme_constraints(model_structure) -> list:
    return(list(model_structure.EnzymeConstraintsInfo.Elements.keys()))


def return_metabolite_constraints(model_structure) -> list:
    return(list(model_structure.MetaboliteConstraintsInfo.Elements.keys()))


def return_model_statistics_information(model_structure) -> list:
    return(model_structure.ModelStatistics.Elements)


def return_general_model_information(model_structure) -> list:
    return(model_structure.GeneralInfo.Elements)


def return_compartment_information(model_structure,compartment: str) -> dict:
    try:
        return(model_structure.CompartmentInfo.Elements[compartment])
    except:
        warnings.warn("{} {} not found".format("Compartment",compartment))


def return_enzyme_information(model_structure,enzyme: str) -> dict:
    try:
        return(model_structure.EnzymeInfo.Elements[enzyme])
    except:
        warnings.warn("{} {} not found".format("Enzyme",enzyme))


def return_protein_information(model_structure,protein: str) -> dict:
    try:
        return(model_structure.ProteinInfo.Elements[protein])
    except:
        warnings.warn("{} {} not found".format("Protein",protein))


def return_reaction_information(model_structure,reaction: str) -> dict:
    try:
        return(model_structure.ReactionInfo.Elements[reaction])
    except:
        warnings.warn("{} {} not found".format("Reaction",reaction))


def return_metabolite_information(model_structure,metabolite: str) -> dict:
    try:
        return(model_structure.MetaboliteInfo.Elements[metabolite])
    except:
        warnings.warn("{} {} not found".format("Metabolite",metabolite))


def return_process_information(model_structure,process: str) -> dict:
    try:
        return(model_structure.ProcessInfo.Elements[process])
    except:
        warnings.warn("{} {} not found".format("Process",process))


def return_macro_molecule_information(model_structure,macro_molecule: str) -> dict:
    try:
        return(model_structure.MacromoleculeInfo.Elements[macro_molecule])
    except:
        warnings.warn("{} {} not found".format("Macromolecule",macro_molecule))


def return_target_information(model_structure,target: str) -> dict:
    try:
        return(model_structure.TargetInfo.Elements[target])
    except:
        warnings.warn("{} {} not found".format("Target",target))


def return_module_information(model_structure,module: str) -> dict:
    try:
        return(model_structure.ModuleInfo.Elements[module])
    except:
        warnings.warn("{} {} not found".format("Module",module))


def return_density_constraint_information(model_structure,density_constraint: str) -> dict:
    try:
        return(model_structure.DensityConstraintsInfo.Elements[density_constraint])
    except:
        warnings.warn("{} {} not found".format("Constraint",density_constraint))


def return_process_constraint_information(model_structure,process_constraint: str) -> dict:
    try:
        return(model_structure.ProcessConstraintsInfo.Elements[process_constraint])
    except:
        warnings.warn("{} {} not found".format("Constraint",process_constraint))


def return_metabolite_constraint_information(model_structure,metabolite_constraint: str) -> dict:
    try:
        return(model_structure.MetaboliteConstraintsInfo.Elements[metabolite_constraint])
    except:
        warnings.warn("{} {} not found".format("Constraint",metabolite_constraint))


def return_enzyme_constraint_information(model_structure,enzyme_constraint: str) -> dict:
    try:
        return(model_structure.EnzymeConstraintsInfo.Elements[enzyme_constraint])
    except:
        warnings.warn("{} {} not found".format("Constraint",enzyme_constraint))


def return_parameter_definition(model,parameter: str) -> dict:
        if parameter in model.parameters.functions._elements_by_id.keys():
            function = model.parameters.functions._elements_by_id[parameter]
            expression = parse_function(function)
            return(expression)
        elif parameter in model.parameters.aggregates._elements_by_id.keys():
            function_id_list = get_function_list_of_aggregate(
                aggregate=model.parameters.aggregates._elements_by_id[parameter])
            expression = parse_aggregate(aggregate=model.parameters.aggregates._elements_by_id[parameter], function_list=[model.parameters.functions._elements_by_id[f_id] for f_id in function_id_list])
            return(expression)
        else:
            return({})


def join_functions_multiplicatively(parsed_function_list: list) -> dict:
    term_list = []
    variable_list = []
    function_type_counts={}
    type_function_dict={}
    function_parameter_list=[]
    multiplicative_term_list=[]
    for function in parsed_function_list:
        function_ID = list(function.keys())[0]
        term_list.append(str('('+function[function_ID]['Equation']+')'))
        variable_list += function[function_ID]['Variables']
        function_parameter_list += function[function_ID]['Function_parameters']
        multiplicative_term_list.append(function_ID)
        if function[function_ID]['Type'] in type_function_dict.keys():
            type_function_dict[function[function_ID]['Type']].append(function)
            function_type_counts[function[function_ID]['Type']]+=1
        else:
            type_function_dict[function[function_ID]['Type']]=[function]
            function_type_counts[function[function_ID]['Type']]=1
    generic_latex_terms=[]
    specific_latex_terms=[]
    for typ in function_type_counts:
        if function_type_counts[typ]==1:
            if typ=='linear':
                generic_latex_terms.append('({})'.format(type_function_dict[typ][0][list(type_function_dict[typ][0].keys())[0]]['Generic_latex'][1:-1]))
                specific_latex_terms.append('({})'.format(type_function_dict[typ][0][list(type_function_dict[typ][0].keys())[0]]['Specific_latex'][1:-1]))
            else:
                generic_latex_terms.append(str(type_function_dict[typ][0][list(type_function_dict[typ][0].keys())[0]]['Generic_latex'][1:-1]))
                specific_latex_terms.append(str(type_function_dict[typ][0][list(type_function_dict[typ][0].keys())[0]]['Specific_latex'][1:-1]))
        else:
            count=0
            for function in type_function_dict[typ]:
                if typ=='linear':
                    generic_latex_terms.append('({})'.format(function[list(function.keys())[0]]['Generic_latex'][1:-1].replace('a','a_{{{}}}'.format(count)).replace('b','b_{{{}}}'.format(count))))
                elif typ=='constant':
                    generic_latex_terms.append(function[list(function.keys())[0]]['Generic_latex'][1:-1].replace('C','C_{{{}}}'.format(count)))
                    specific_latex_terms.append(str(function[list(function.keys())[0]]['Specific_latex'][1:-1]))
                elif typ=='exponential':
                    generic_latex_terms.append(function[list(function.keys())[0]]['Generic_latex'][1:-1].replace('_/_','_/__{{{}}}'.format(count)))
                    specific_latex_terms.append(str(function[list(function.keys())[0]]['Specific_latex'][1:-1]))
                elif typ=='michaelisMenten':
                    generic_latex_terms.append(function[list(function.keys())[0]]['Generic_latex'][1:-1].replace('V','V_{{{}}}'.format(count)).replace('K','K_{{{}}}'.format(count)))
                    specific_latex_terms.append(str(function[list(function.keys())[0]]['Specific_latex'][1:-1]))
                count+=1
    return({'Type': 'Aggregate',
            'Equation': '*'.join(term_list),
            'Variables': list(set(variable_list)),
            'Generic_latex':'${}$'.format(' _/_times '.join(generic_latex_terms)),
            'Specific_latex':'${}$'.format(' _/_times '.join(specific_latex_terms)),
            'Function_parameters':function_parameter_list,
            'Multiplicative Terms':multiplicative_term_list})


def parse_aggregate(aggregate, function_list: list) -> dict:
    aggregate_ID = aggregate.id
    if aggregate.type == 'multiplication':
        parsed_function_list = [parse_function(function) for function in function_list]
        result = {aggregate_ID: join_functions_multiplicatively(parsed_function_list=parsed_function_list)}
        return(result)
    else:
        return({aggregate_ID: {'Type': 'Aggregate',
                               'Equation': '',
                               'Variables': [],
                               'Generic_latex': '',
                               'Specific_latex': '',
                               'Multiplicative Terms': []}})


def parse_function(function) -> dict:
    independent_variable = function.variable
    if independent_variable=='growth_rate':
        latex_independent_variable='_/_mu'
    elif independent_variable.endswith('_e'):
        if independent_variable.startswith('M_'):
            latex_independent_variable='c_{{{}}}'.format(independent_variable[2:-2])
        else:
            latex_independent_variable='c_{{{}}}'.format(independent_variable[:-2])
    else:
        if independent_variable.startswith('M_'):
            latex_independent_variable='c_{{{}}}'.format(independent_variable[2:])
        else:
            latex_independent_variable='c_{{{}}}'.format(independent_variable)

    function_ID = function.id
    if function.type == 'constant':
        eq = make_paramter_function_specific(
            function_ID=function_ID, parameter='CONSTANT', return_normal=True)
        latex_string_generic = '$C$'
        latex_string_specific = '${}$'.format(get_function_parameter_value_from_model(function=function, parameter_ID='CONSTANT'))
        function_parameter_values = {'CONSTANT': get_function_parameter_value_from_model(function=function, parameter_ID='CONSTANT')}

    elif function.type == 'exponential':
        eq = 'e**({}*{})'.format(make_paramter_function_specific(function_ID=function_ID,
                                                                 parameter='RATE', return_normal=True), str(independent_variable))
        latex_string_generic = '$e^{{_/_lambda {}}}$'.format(latex_independent_variable)
        latex_string_specific = '$e^{{{} {}}}$'.format(get_function_parameter_value_from_model(function=function, parameter_ID='RATE'),latex_independent_variable)
        function_parameter_values = {'e':numpy.e,'RATE': get_function_parameter_value_from_model(function=function, parameter_ID='RATE')}

    elif function.type == 'linear':
        eq = str('{}+{}*{}'.format(make_paramter_function_specific(function_ID=function_ID, parameter='LINEAR_CONSTANT', return_normal=True),
                                   make_paramter_function_specific(function_ID=function_ID, parameter='LINEAR_COEF', return_normal=True), str(independent_variable)))
        latex_string_generic ='$a {} + b$'.format(latex_independent_variable)
        latex_string_specific ='${} {} + {}$'.format(get_function_parameter_value_from_model(function=function, parameter_ID='LINEAR_COEF'),latex_independent_variable,get_function_parameter_value_from_model(function=function, parameter_ID='LINEAR_CONSTANT'))
        function_parameter_values = {'LINEAR_CONSTANT': get_function_parameter_value_from_model(function=function, parameter_ID='LINEAR_CONSTANT'),
                                     'LINEAR_COEF': get_function_parameter_value_from_model(function=function, parameter_ID='LINEAR_COEF'),
                                     'X_MIN': get_function_parameter_value_from_model(function=function, parameter_ID='X_MIN'),
                                     'X_MAX': get_function_parameter_value_from_model(function=function, parameter_ID='X_MAX'),
                                     'Y_MIN': get_function_parameter_value_from_model(function=function, parameter_ID='Y_MIN'),
                                     'Y_MAX': get_function_parameter_value_from_model(function=function, parameter_ID='Y_MAX'), }

    elif function.type == 'michaelisMenten':
        eq = str('{}*{}/({}+{})'.format(make_paramter_function_specific(function_ID=function_ID, parameter='kmax', return_normal=True),
                                        str(independent_variable), str(independent_variable), make_paramter_function_specific(function_ID=function_ID, parameter='Km', return_normal=True)))
        latex_string_generic ='$_/_frac{{V{}}}{{{}+K}}$'.format(latex_independent_variable,latex_independent_variable)
        latex_string_specific ='$_/_frac{{{}{}}}{{{}+{}}}$'.format(get_function_parameter_value_from_model(function=function, parameter_ID='kmax'),latex_independent_variable,latex_independent_variable,get_function_parameter_value_from_model(function=function, parameter_ID='Km'))
        function_parameter_values = {'kmax': get_function_parameter_value_from_model(function=function, parameter_ID='kmax'),
                                     'Km': get_function_parameter_value_from_model(function=function, parameter_ID='Km'),
                                     'Y_MIN': get_function_parameter_value_from_model(function=function, parameter_ID='Y_MIN')}

    return({function_ID: {'Type': function.type,
                          'Generic_latex':latex_string_generic,
                          'Specific_latex':latex_string_specific,
                          'Equation': eq,
                          'Variables': [str(independent_variable)],
                          'Function_parameters': function_parameter_values,
                          'Multiplicative Terms':[function_ID]}})


def make_paramter_function_specific(function_ID, parameter, return_normal=False) -> str:
    if return_normal:
        return(str(parameter))
    else:
        return(str('{}__parameter__{}'.format(function_ID, parameter)))


def get_function_list_of_aggregate(aggregate) -> list:
    return([agg.function for agg in aggregate.function_references._elements])


def get_function_parameter_value_from_model(function, parameter_ID) -> float:
    if parameter_ID in list(function.parameters._elements_by_id.keys()):
        return(function.parameters._elements_by_id[parameter_ID].value)
    else:
        return(None)


def evaluate_current_function_value(expression_dictionary: dict , growth_rate: float , medium: dict) -> float:
    eq=list(expression_dictionary.values())[0]['Equation']
    function_params=list(expression_dictionary.values())[0]['Function_parameters']
    for v in list(expression_dictionary.values())[0]['Variables']:
        if v == 'growth_rate':
            function_params[v] = growth_rate
            if 'X_MAX' in function_params:
                if growth_rate>function_params['X_MAX']:
                    function_params[v] = function_params['X_MAX']
            if 'X_MIN' in function_params:
                if growth_rate<function_params['X_MIN']:
                    function_params[v] = function_params['X_MIN']
        elif v in medium.keys():
            function_params[v] = medium[v]
            if 'X_MAX' in function_params:
                if medium[v]>function_params['X_MAX']:
                    function_params[v] = function_params['X_MAX']
            if 'X_MIN' in function_params:
                if medium[v]<function_params['X_MIN']:
                    function_params[v] = function_params['X_MIN']
        elif v.endswith('_e'):
            if v[:-2] in medium.keys():
                function_params[v] = medium[v[:-2]]
                if 'X_MAX' in function_params:
                    if medium[v[:-2]]>function_params['X_MAX']:
                        function_params[v] = function_params['X_MAX']
                if 'X_MIN' in function_params:
                    if medium[v[:-2]]<function_params['X_MIN']:
                        function_params[v] = function_params['X_MIN']
        else:
            return(numpy.nan)
    result=eval(eq,function_params)
    if 'Y_MAX' in function_params:
        if function_params['Y_MAX'] is not None:
            if result>function_params['Y_MAX']:
                return(function_params['Y_MAX'])
            else:
                return(result)
        else:
            return(result)
    elif 'Y_MIN' in function_params:
        if function_params['Y_MIN'] is not None:
            if result<function_params['Y_MIN']:
                return(function_params['Y_MIN'])
            else:
                return(result)
        else:
            return(result)
    else:
        return(result)


def get_medium_dependent_coefficients_in_lhs(ModelStructure):
    medium_dependent_constraints=[]
    for met in ModelStructure.MediumDependencies.keys():
        medium_dependent_constraints+=ModelStructure.MediumDependencies[met]
    out=[]
    for i in list(set(medium_dependent_constraints)):
        if i in ModelStructure.EnzymeConstraintsInfo.Elements.keys():
            out.append((i,ModelStructure.EnzymeConstraintsInfo.Elements[i]['AssociatedEnzyme']))
        elif i in ModelStructure.ProcessConstraintsInfo.Elements.keys():
            out.append((i,str(ModelStructure.ProcessConstraintsInfo.Elements[i]['AssociatedProcessID']+"_machinery")))
    return(out)

def record_proto_proteome(RBA_Session, run, Proteinlevels):
    out = []
    for i in list(RBA_Session.ModelStructure.ProteinGeneMatrix['ProtoProteins']):
        row_ind = list(RBA_Session.ModelStructure.ProteinGeneMatrix['ProtoProteins']).index(i)
        nonZero = list(numpy.nonzero(RBA_Session.ModelStructure.ProteinGeneMatrix['Matrix'][row_ind, :])[0])
        level = 0
        for j in nonZero:
            id = RBA_Session.ModelStructure.ProteinGeneMatrix['Proteins'][j]
            level += Proteinlevels.loc[id, run]
        out.append(level)
    return(out)


def record_proteome(RBA_Session, run):

    EnzDF = pandas.DataFrame(index=RBA_Session.Problem.Enzymes)
    PrcDF = pandas.DataFrame(index=RBA_Session.Problem.Processes)
    EnzDF[run] = [RBA_Session.Problem.SolutionValues[i]for i in RBA_Session.Problem.Enzymes]
    PrcDF[run] = [RBA_Session.Problem.SolutionValues[i]for i in RBA_Session.Problem.Processes]

    ProteinProteinMatrix = numpy.array(
        RBA_Session.ModelStructure.ProteinMatrix['Matrix']).astype(numpy.float64)
    C = RBA_Session.ModelStructure.ProteinMatrix['Consumers']
    Consumers = []
    for i in C:
        if i.startswith('P_'):
            Consumers.append(str(i))
        if not i.startswith('P_'):
            Consumers.append(i)
    Proteins = RBA_Session.ModelStructure.ProteinMatrix['Proteins']
    DF = pandas.concat([EnzDF, PrcDF], axis=0)
    ProteinLevels = pandas.DataFrame(index=Proteins)
    vector = numpy.nan_to_num(DF[run].reindex(Consumers))
    Level = ProteinProteinMatrix.dot(vector)
    ProteinLevels[run] = Level
    addedProts = [col for col in RBA_Session.Problem.LP.col_names if col.startswith('TotalLevel_')]
    if len(addedProts) > 0:
        for p in addedProts:
            protID = p.split('TotalLevel_')[1]
            ProteinLevels[run].loc[protID] = RBA_Session.Problem.SolutionValues[p]
    return(list(ProteinLevels[run]))


def map_iso_reactions(RBA_Session):
    if hasattr(RBA_Session, 'Results'):
        out = pandas.DataFrame()
        for run in list(RBA_Session.Results['Reactions'].columns):
            rf = dict(zip(list(RBA_Session.Results['Reactions'].index), list(
                RBA_Session.Results['Reactions'][run])))
            rf = {k: v for k, v in rf.items() if v != 0.}
            rf_merged = collections.defaultdict(float)
            for reac_id, flux_val in rf.items():
                if "duplicate" in reac_id:
                    last_idx = reac_id.index('duplicate') - 1
                    rf_merged[reac_id[:last_idx]] += flux_val
                else:
                    rf_merged[reac_id] += flux_val
            if not list(out):
                out[run] = list(rf_merged.values())
                out.index = list(rf_merged.keys())
            else:
                runDF = pandas.DataFrame(list(rf_merged.values()),
                                         index=list(rf_merged.keys()), columns=[run])
                runDF = runDF.reindex(list(set(list(out.index)).union(
                    set(list(rf_merged.keys())))), fill_value=0)
                out = out.reindex(list(set(list(out.index)).union(
                    set(list(rf_merged.keys())))), fill_value=0)
                out = out.join(runDF, how='outer')
        return(out)


def build_exchange_map(RBA_Session):
    """
    Returns a map of all metabolites, the corresponding transport-reactions and stoichiometires;
    exchanged with the medium.
    {Metabolite1 : {ExchangeReaction1 : stoch-coefficient1 , ExchangeReaction2 : stoch-coefficient2},
    {Metabolite2 : {ExchangeReaction1 : stoch-coefficient1 , ExchangeReaction2 : stoch-coefficient2}}

    Metabolite1 - ... MetaboliteN : All metabolite-species in the medium (see medium.tsv file)
    ExchangeReaction1 - ... ExchangeReactionN : All metabolic reactions, which exchange the respective metabolite with the medium.
    stoch-coefficient : Stochiometric coefficient with which the respective metabolite is exchanged by the corresponding reaction.
    (Negative when reaction transports metabolite out of the cell; and positive when inside the cell.)

    Parameters
    ----------
    RBA_Session : rbatools.NewControler.RBA_newControler

    Returns
    -------
    Dict.
    """
    BoundaryMetabolites = [i for i in list(RBA_Session.ModelStructure.MetaboliteInfo.Elements.keys(
    )) if RBA_Session.ModelStructure.MetaboliteInfo.Elements[i]['boundary']]
    ExchangeMap = {}
    for bM in BoundaryMetabolites:
        for rxn in RBA_Session.ModelStructure.MetaboliteInfo.Elements[bM]['ReactionsInvolvedWith']:
            Reactants = list(
                RBA_Session.ModelStructure.ReactionInfo.Elements[rxn]['Reactants'].keys())
            Products = list(RBA_Session.ModelStructure.ReactionInfo.Elements[rxn]['Products'].keys())
            if len(list(set(list(Reactants+Products)))) > 1:
                for met in list(set(list(Reactants+Products))):
                    # if met != bM:
                    if met == bM:
                        MediumSpecies = find_exchange_metabolite_in_medium(metabolite=met, Medium=RBA_Session.Medium)
                        if met in Reactants:
                            stochCoeff = - \
                                RBA_Session.ModelStructure.ReactionInfo.Elements[rxn]['Reactants'][met]
                        elif met in Products:
                            stochCoeff = RBA_Session.ModelStructure.ReactionInfo.Elements[rxn]['Products'][met]
                        if MediumSpecies in list(ExchangeMap.keys()):
                            ExchangeMap[MediumSpecies].update({rxn: stochCoeff})
                        else:
                            ExchangeMap[MediumSpecies] = {rxn: stochCoeff}
    return(ExchangeMap)


def find_exchange_metabolite_in_medium(metabolite, Medium):
    """
    Returns the most likely species in the Medium, for any Metabolic species.
    Parameters
    ----------
    metabolite : str
    Medium : dict
    -------
    Most likely ID as str
    """
    if metabolite.endswith('_e'):
        out = difflib.get_close_matches('_e'.join(metabolite.split('_e')[:-1]), Medium, 1)
    else:
        out = difflib.get_close_matches(metabolite, Medium, 1)
    if len(out) > 0:
        return(out[0])
    else:
        return('')


def check_solution_feasibility(Value, RBA_Session):
    if RBA_Session.Problem.Solved:
        return(Value)
    else:
        return(numpy.nan)

def convert_cplex_matrix_to_sparse(inputStructure):
    Ma = inputStructure.cplexLP.linear_constraints.get_rows()
    Anew = numpy.zeros((inputStructure.cplexLP.linear_constraints.get_num(),
                        inputStructure.cplexLP.variables.get_num()))
    rowIndex = 0
    for m in Ma:
        Anew[rowIndex, m.ind] = m.val
        rowIndex += 1
    return(scipy.sparse.coo_matrix(Anew))

def check_for_attributes(obj, attr):
    for i in attr:
        x = getattr(obj, i, None)
        if x is None:
            return(True)
    return(True)
