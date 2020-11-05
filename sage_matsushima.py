from scipy import stats
from scipy.stats import t as t_dist
from scipy.stats import chi2

import math
import statistics
from abtesting_test import *

# # You can comment out these lines! They are just here to help follow along to the tutorial.
# print(t_dist.cdf(-2, 20))  # should print .02963
# # positive t-score (bad), should print .97036 (= 1 - .2963)
# print(t_dist.cdf(2, 20))

# print(chi2.cdf(23.6, 12))  # prints 0.976
# print(1 - chi2.cdf(23.6, 12))  # prints 1 - 0.976 = 0.023 (yay!)

# TODO: Fill in the following functions! Be sure to delete "pass" when you want to use/run a function!
# NOTE: You should not be using any outside libraries or functions other than the simple operators (+, **, etc)
# and the specifically mentioned functions (i.e. round, cdf functions...)


def slice_2D(list_2D, start_row, end_row, start_col, end_col):
    '''
    Splices a the 2D list via start_row:end_row and start_col:end_col
    :param list: list of list of numbers
    :param nums: start_row, end_row, start_col, end_col
    :return: the spliced 2D list (ending indices are exclsive)
    '''
    to_append = []
    for l in range(start_row, end_row):
        to_append.append(list_2D[l][start_col:end_col])

    return to_append


def get_avg(nums):
    '''
    Helper function for calculating the average of a sample.
    :param nums: list of numbers
    :return: average of list
    '''
    sum = 0
    for numbers in nums:
        sum += numbers

    return sum/len(nums)


def get_stdev(nums):
    '''
    Helper function for calculating the standard deviation of a sample.
    :param nums: list of numbers
    :return: standard deviation of list
    '''
    n = len(nums)

    if n <= 1:
        return 0.0

    mean = get_avg(nums)
    sd = 0

    for numbers in nums:
        sd += (float(numbers) - mean)**2
    sd = math.sqrt(sd / float(n-1))

    return sd


def get_standard_error(a, b):
    '''
    Helper function for calculating the standard error, given two samples.
    :param a: list of numbers
    :param b: list of numbers
    :return: standard error of a and b (see studio 6 guide for this equation!)
    '''
    return math.sqrt(((get_stdev(a)**2)/len(a)) + ((get_stdev(b)**2)/len(b)))


def get_2_sample_df(a, b):
    '''
    Calculates the combined degrees of freedom between two samples.
    :param a: list of numbers
    :param b: list of numbers
    :return: integer representing the degrees of freedom between a and b (see studio 6 guide for this equation!)
    HINT: you can use Math.round() to help you round!
    '''
    se = get_standard_error(a, b)
    stda = get_stdev(a)**2/len(a)
    stdb = get_stdev(b)**2/len(b)

    return round((se**4)/(((stda**2)/(len(a) - 1)) + ((stdb**2)/(len(b) - 1))))


def get_t_score(a, b):
    '''
    Calculates the t-score, given two samples.
    :param a: list of numbers
    :param b: list of numbers
    :return: number representing the t-score given lists a and b (see studio 6 guide for this equation!)
    '''
    return ((get_avg(a) - get_avg(b))/get_standard_error(a, b))


def perform_2_sample_t_test(a, b):
    '''
    ** DO NOT CHANGE THE NAME OF THIS FUNCTION!! ** (this will mess with our autograder)
    Calculates a p-value by performing a 2-sample t-test, given two lists of numbers.
    :param a: list of numbers
    :param b: list of numbers
    :return: calculated p-value
    HINT: the t_dist.cdf() function might come in handy!
    '''
    return t_dist.cdf(get_t_score(a, b), get_2_sample_df(a, b))


# [OPTIONAL] Some helper functions that might be helpful in get_expected_grid().
def row_sum(observed_grid, ele_row):
    sum = 0
    new = slice_2D(observed_grid, ele_row, ele_row+1, 0, len(observed_grid[0]))
    for numbers in range(len(new[0])):
        sum += new[0][numbers]
    return sum


def col_sum(observed_grid, ele_col):
    sum = 0
    new = slice_2D(observed_grid, 0, len(observed_grid), ele_col, ele_col+1)
    for numbers in range(len(new[0])+1):
        sum += new[numbers][0]
    return sum


def total_sum(observed_grid):
    sum = 0
    for row in range(len(observed_grid)):
        for col in range(len(observed_grid[row])):
            sum += observed_grid[row][col]
    return sum


def calculate_expected(row_sum, col_sum, tot_sum):
    return ((row_sum*col_sum)/tot_sum)


def get_expected_grid(observed_grid):
    '''
    Calculates the expected counts, given the observed counts.
    ** DO NOT modify the parameter, observed_grid. **
    :param observed_grid: 2D list of observed counts
    :return: 2D list of expected counts
    HINT: To clean up this calculation, consider filling in the optional helper functions below!
    '''
    sum = total_sum(observed_grid)
    expected = []

    for row in range(0, len(observed_grid)):
        x_row = []
        for col in range(0, len(observed_grid[row])):
            rsum = row_sum(observed_grid, row)
            csum = col_sum(observed_grid, col)
            exval = calculate_expected(rsum, csum, sum)
            x_row.append(exval)
        expected.append(x_row)

    return expected


def df_chi2(observed_grid):
    '''
    Calculates the degrees of freedom of the expected counts.
    :param observed_grid: 2D list of observed counts
    :return: degrees of freedom of expected counts (see studio 6 guide for this equation!)
    '''
    rows = len(observed_grid) - 1
    col = len(observed_grid[0]) - 1
    return (rows*col)


def chi2_value(observed_grid):
    '''
    Calculates the chi^2 value of the expected counts.
    :param observed_grid: 2D list of observed counts
    :return: associated chi^2 value of expected counts (see studio 6 guide for this equation!)
    '''
    chi2 = 0
    expected = get_expected_grid(observed_grid)

    for row in range(0, len(expected)):
        for col in range(0, len(expected[row])):
            chi2 += ((observed_grid[row][col] -
                      expected[row][col])**2)/(expected[row][col])
    return chi2


def perform_chi2_homogeneity_test(observed_grid):
    '''
    ** DO NOT CHANGE THE NAME OF THIS FUNCTION!! ** (this will mess with our autograder)
    Calculates the p-value by performing a chi^2 test, given a list of observed counts
    :param observed_grid: 2D list of observed counts
    :return: calculated p-value
    HINT: the chi2.cdf() function might come in handy!
    '''
    return 1 - chi2.cdf(chi2_value(observed_grid), df_chi2(observed_grid))

# These commented out lines are for testing your main functions.
# Please uncomment them when finished with your implementation and confirm you get the same values :)


def data_to_num_list(s):
    '''
      Takes a copy and pasted row/col from a spreadsheet and produces a usable list of nums.
      This will be useful when you need to run your tests on your cleaned log data!
      :param str: string holding data
      :return: the spliced list of numbers
      '''
    return list(map(float, s.split()))


# chi2_test 3:
ca_list = data_to_num_list(ac)
cb_list = data_to_num_list(bc)
cog = [ca_list, cb_list]
print(chi2_value(cog))  # 0.05099067599067608
print(perform_chi2_homogeneity_test(cog))  # 0.8213483036962936

# T TEST:
ta_list = data_to_num_list(a)
tb_list = data_to_num_list(b)
print(get_t_score(ta_list, tb_list))  # this should be -0.643263379330903
# this should be 0.2666168298376804
print(perform_2_sample_t_test(ta_list, tb_list))
