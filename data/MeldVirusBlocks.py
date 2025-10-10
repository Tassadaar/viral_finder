import os
import pandas as pd
import gffutils
import sqlite3
import PandasSettings
import pprint

PandasSettings.Set_Panda_Display()

search_radius = 20000

def id_func(f) -> str:
    '''
    If merge_strategy='create_unique' fails,
    use this function to create unique id's
    as a fallback strategy.

    Particularly useful for Liftoff-generated
    GFF3's, which like gffutils also appends
    _1 to identical feature IDs
    '''
    if f.featuretype in ['gene', 'mRNA']:
        return f.attributes['ID'][0]
    elif f.featuretype in ['exon', 'CDS']:
        return '-'.join( [f.attributes['ID'][0], f.seqid, str(f.start), str(f.end)] )

def get_gff(gff):
    try:
        gff3_db = gffutils.create_db(gff, dbfn=':memory:', merge_strategy='create_unique')
    except sqlite3.IntegrityError:
        print("Default loading strategy 'create_unique' failed,",
              "changing to custom id_spec function strategy", )
        gff3_db = gffutils.create_db(gff, dbfn=':memory:', id_spec=id_func)
    return gff3_db

os.chdir("J:/Assembly_and_Genes_sept2025/GenePredictionsSept2025/ST7C/Final_Gff3s/")

db = get_gff("ST7C_gene_predictions_viraltag.gff3")

#convert to table
gff_to_table =[]
for f in db.features_of_type(featuretype='gene'):
    gff_row = {"seqid":"", "genestart":"", "genestop":"", "ID":"", "Protein type":"", "class":""}

    gff_row['seqid'] = f.seqid
    gff_row['genestart'] = f.start
    gff_row['genestop'] = f.end
    gff_row['ID'] = f.attributes['ID'][0]
    # type_of_protein = f.attributes['Name'][0]
    if 'Name' in f.attributes:
        classtype = 'viral'
        type_of_protein = f.attributes['Name'][0].split("MELDVirus_")[-1].replace(" ", "_").replace(",", "")
    else:
        classtype = 'Non_viral'
        type_of_protein = "Non_viral"

    gff_row['Protein type'] = type_of_protein
    gff_row['class'] = classtype
    gff_to_table.append(gff_row)

df = pd.DataFrame(gff_to_table)
# df = df.set_index('contig')
pd.set_option('display.max_rows', None)
print(df)
exit()
contigs_to_look_at = sorted(list(set(df['seqid'].to_list())))
# print(contigs_to_look_at)
# list_of_viral_blocks=[]

for contig in contigs_to_look_at:
    sub_df = df[df['seqid'] == contig].copy()

    # print(sub_df)

    sub_df['block'] = sub_df['class'] != sub_df['class'].shift()
    sub_df['run_id'] = sub_df['block'].cumsum()
    # print(sub_df)
    max_gap = 2

    runs = sub_df.groupby('run_id').agg(
        start=('run_id', 'idxmin'),
        end=('run_id', 'idxmax'),
        viral_block = ('class', 'first')
    ).reset_index(drop=True)

    merged = []
    current = runs.iloc[0].to_dict()
    print(runs)
    for _, r in runs.iloc[1:].iterrows():
        if (
            r['viral_block'] == current['viral_block'] and
            r['start'] - current['end'] - 1 <= max_gap
        ):
            current['end'] = r['end']
        else:
            merged.append(current)
            current = r.to_dict()
    merged.append(current)
    # print(contig)
    pprint.pprint(merged)

    groups = [sub_df.loc[m['start']:m['end']] for m in merged]
    groups = [x for x in groups if len(x)> 1]
    # print(groups)


# for contig in contigs_to_look_at:
#     sub_df = df.loc[contig]
#     sub_df = sub_df.reset_index()
#     entries_examined = []
#     for i, row in sub_df.iterrows():
#         print(row)
#         start_value = row['start']
#         internal_non_virals = 0
#         block_of_genes = []
#         if row['Protein type'] != 'Non_viral' and not row['ID'] in entries_examined:
#             entries_examined.append(row['ID'])
#             block_of_genes.append(row['ID'])
#             continue_bool = True
#             while continue_bool:
#                 print(i+1)
#                 next_row = sub_df.iloc[i+1]
#                 if next_row['Protein type'] != 'Non_viral':
#                     i = i + 1
#                     entries_examined.append(next_row['ID'])
#                     block_of_genes.append(next_row['ID'])
#                 elif next_row['Protein type'] == 'Non_viral':
#                     next_next_row = sub_df.iloc[i+2]
#                     if next_next_row['Protein type'] != 'Non_viral':
#                         entries_examined.append(next_row['ID'])
#                         i = i +1
#                     else:
#                         continue_bool = False
#         list_of_viral_blocks.append(block_of_genes)


# print(list_of_viral_blocks)



