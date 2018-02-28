import functools
import itertools
import sympy


def set2ls(s, stop):
    converted = sorted(list(s)) + [stop]
    return [y - x for x, y in zip(converted, converted[1:])]


def setrange(chroma=1, stop=10):
    return set(range(0, stop, chroma))


def superunion(*sets):
    r = set([])
    for s in sets:
        r = r | s
    return r


def synchronize(*generators, convert_set2ls=True):
    """
    basic synchronisation:
        * the result is symmetric
        * the results sum is: generator * generator * ... * generator
    """
    cp = functools.reduce(lambda x, y: x * y, generators)  # common_product
    gen = (setrange(generator, cp) for generator in generators)
    union = superunion(*gen)
    if convert_set2ls is True:
        return set2ls(union, cp)
    else:
        return union


def synchronize_complementary(*generators):
    """
    the results length is sum(generators) - 2
    """
    cp = functools.reduce(lambda x, y: x * y, generators)  # common_product
    # complementary factors
    com_fac = (cp // generator for generator in generators)
    gen = (setrange(com_generator, cp) for com_generator in com_fac)
    return set2ls(superunion(*gen), cp)


def fractionize(*generators):
    def mk_fractionized_gen(generator, major, stop):
        new_gen = setrange(generator, generator * major)
        gen = (set([el + (major * counter) for el in new_gen])
               for counter in range(major - generator + 1))
        return superunion(*gen)
    major = max(generators)
    stop = pow(major, 2)  # since complementary_factor of "a" is "a" -> a * a
    gen = (mk_fractionized_gen(generator, major, stop)
           for generator in filter(lambda x: x < major, generators))
    return set2ls(superunion(*gen, setrange(major, stop)), stop)


def mk_complementary_factors(*generators):
    """
    makes complementary_factors of
    input generators
    """
    cp = functools.reduce(lambda x, y: x * y, generators)  # common_product
    return list(cp // generator for generator in generators)


def cyclic(iterable):
    am = len(iterable)
    return tuple(tuple(iterable[i - j] for i in range(am)) for j in range(am))


def permute_general(iterable, lv, function):
    """
    general higher level permutation function
    """
    for counter in range(lv):
        iterable = function(iterable)
        if counter > 0:
            iterable = (functools.reduce(lambda x, y: x + y, element)
                        for element in iterable)
    return iterable


def permute(iterable, lv=1):
    return permute_general(iterable, lv, itertools.permutations)


def permute_cyclic(iterable, lv=1):
    return permute_general(iterable, lv, cyclic)


def distributive_power(*args, power=2):
    class pairedint(int):
        mate = None

    def seperate_coeff(big_term):
        coefficients = []
        terms = []
        for term in big_term:
            if term.func == sympy.Pow:
                coefficients.append(1)
                mul2 = [term.args[0]] * term.args[1]
            else:
                mul2 = []
                for arg in term.args:
                    if arg.is_Integer is True:
                        coefficients.append(arg)
                    else:
                        mul2.append(arg)
            terms.append(sympy.Mul(*mul2))
        return (terms, coefficients)

    def term2args(*terms):
        res = []
        for t in terms:
            args = []
            if t.func == sympy.Pow:
                args += [t.args[0]] * t.args[1]
            else:
                for arg in t.args:
                    if arg.func == sympy.Pow:
                        args += [arg.args[0]] * arg.args[1]
                    else:
                        args += [arg]
            res.append(args)
        return res

    def count_coefficients(sep_term, coefficients):
        am = len(coefficients)
        res = [0] * am
        for arg in sep_term:
            for counter, coeff in enumerate(coefficients):
                if arg == coeff:
                    res[counter] += 1
                    break
        return res

    def distribute_products(terms, term_coefficents,
                            integers, coefficients, power):
        coefficient_groups = [[] for i in range(len(coefficients))]
        for t, t_coeff, i in zip(terms, term_coefficents, integers):
            w = i / power
            for counter, coeff_amount in enumerate(t_coeff):
                coefficient_groups[counter] += [t] * int(w * coeff_amount)
        return coefficient_groups

    def calc_terms(terms, coefficients, args):
        subls = [(coeff, arg) for coeff, arg in zip(coefficients, args)]
        return tuple(term.subs(subls) for term in terms)

    def sort_coefficent_groups(groups, terms, term_coefficents, coefficients):
        def value_coefficents(term_coefficents, coefficients):
            res = []
            for t in term_coefficents:
                v = 0
                for counter, coeff in enumerate(t):
                    v += (counter + 1) * coeff
                res.append(v)
            return res
        term_values = value_coefficents(term_coefficents, coefficients)
        new_groups = []
        for group in groups:
            gr_values = []
            for t in group:
                for counter, t_comp in enumerate(terms):
                    if t == t_comp:
                        val = pairedint(term_values[counter])
                        val.mate = t
                        gr_values.append(val)
                        break
            new_groups.append([val.mate for val in sorted(gr_values)])
        return new_groups

    """ugly implementation of schillingers distributive_power-algorithm"""

    if power > 1:
        var = " + ".join("x{0}".format(counter)
                         for counter in range(len(args)))
        var_term = sympy.sympify(var)
        coefficients = var_term.args
        term = sympy.Pow(var_term, power)
        expanded = term.expand()

        # seperate single terms from their coefficients
        terms, integers = seperate_coeff(expanded.args)

        # calculate single terms
        terms_calculated = calc_terms(terms, coefficients, args)

        # seperate terms to elements
        seperated_term = term2args(*terms)

        # count coefficients of every term
        term_coefficents = [count_coefficients(
            t, coefficients) for t in seperated_term]

        # distribute different terms to coefficient-groups
        coefficient_groups = distribute_products(
            terms, term_coefficents, integers, coefficients, power)

        # sort coefficient_groups
        sorted_coefficent_groups = sort_coefficent_groups(
                coefficient_groups, terms, term_coefficents, coefficients)

        # replace terms with their results
        res = []
        for group in sorted_coefficent_groups:
            for t in group:
                for counter, comp_t in enumerate(terms):
                    if comp_t == t:
                        res.append(terms_calculated[counter])
                        break

        return res

    else:
        return list(args)


def adjust(rh, relation):
    factor = sum(relation) / sum(rh)
    return [el * factor for el in rh]
