import gspread
from oauth2client.service_account import ServiceAccountCredentials

scope = ['https://spreadsheets.google.com/feeds',
         'https://www.googleapis.com/auth/drive']

credentials = ServiceAccountCredentials.from_json_keyfile_name('../VCN-SBEM-Data/client_secret.json', scope)

gc = gspread.authorize(credentials)


allss = gc.open("Dendrite Quality and Surface Areas_comparisons_6_March_2020")
print(dir(allss))
print(allss.worksheets())
wks = {}
for i,  w in enumerate(allss.worksheets()):
    print(w.title)
    wks[w.title] = allss.get_worksheet(i)

print(wks)

# wks1 = gc.open("Dendrite Quality and Surface Areas_comparisons_6_March_2020").sheet1
# wks2 = gc.open("Dendrite Quality and Surface Areas_comparisons_6_March_2020").sheet2
# wks3 = gc.open("Dendrite Quality and Surface Areas_comparisons_6_March_2020").sheet3
#
# listofhashes = wks1.get_all_records()
# print(listofhashes)

