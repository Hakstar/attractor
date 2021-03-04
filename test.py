import itertools
test_list=[1,2,3,4]
a = list(itertools.combinations(test_list,2))
print(a[0][1])