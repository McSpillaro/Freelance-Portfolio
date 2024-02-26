start_year = 2018 # year to start regression from
end_year = start_year + 16 # creating a regression of 17 years

years_months_list_strings = [] # For x-labels on graph to be formatted correctly
tmp_month = 1 # month number for string
tmp_year = start_year # year number for string

while tmp_year != end_year:
    years_months_list_strings.append(f'{tmp_year}-{tmp_month}')
    if tmp_month == 12:
        tmp_year += 1 # goes to next year after last month
        tmp_month = 1 # resets month to 1st one