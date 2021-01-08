with open("nname_field.dat", "r") as f:
    all_fields = f.readlines()

with open("nname_field_avail.dat", "r") as f:
    avai_fields = f.readlines()

with open("nname_field_raw_avail.dat", "r") as f:
    avai_raw_fields = f.readlines()

print("Totaly: %d fields"%len(all_fields))
print("Finaly: %d (%d raw) fields"%(len(avai_fields),len(avai_raw_fields)))

not_avail = []
not_avail_raw = []
for nm in all_fields:
    if nm not in avai_fields:
        not_avail.append(nm)

    if nm not in avai_raw_fields:
        not_avail_raw.append(nm)
print("Not in final cata: ", not_avail)
print("Not in final raw cata: ", not_avail_raw)