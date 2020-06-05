import pandas as pd
import argparse
import pickle
import os


def filter_biobank_af(x):
    # print(x['TaiwanBioBank'])
    if x == '.':
        return 0
    else:
        biobank_af = x.split('|')[2]
        # print(biobank_af)
        AF = biobank_af[3:]
        # print(AF)
        return AF


def filter(source_df):
    # ------- Actionable ---------
    actionable_df = source_df[~source_df['oncoKB_annotation'].isin(['.']) |
                              ~source_df['CGI_annotation'].isin(['.']) |
                              ~source_df['CIVIC_annotation'].isin(['.'])]
    print("------- Actionable variants ---------")
    # print(actionable_df)

    # ------- Filter ---------
    # AF <= 0.01
    tmp_df = source_df[source_df['AF'] <= 0.01]
    # Biobank AF <= 0.01
    tmp_df['biobank_af'] = tmp_df['TaiwanBioBank'].apply(filter_biobank_af)
    # print(df['biobank_af'])
    df = tmp_df[tmp_df['biobank_af'].astype(float) <= 0.01]
    df = df.drop(columns=['biobank_af'])
    print(df)

    # Exonic
    func_list = ['exonic', 'splicing', 'exonic;splicing']
    df = df[df['Func.refGene'].isin(func_list)]

    # Nonsynonymous
    filter_df = df[~df['ExonicFunc.refGene'].isin(['synonymous SNV'])]

    # ----- Heredity -------
    tmp_heredity_df = filter_df[(filter_df['CLNREVSTAT'].isin(['reviewed_by_expert_panel']) &
                                 filter_df['CLNSIG'].isin(['Pathogenic', 'Likely_pathogenic'])) |
                                (filter_df['LOVD_all_clinical'].str.contains('pathogenic') &
                                 (~filter_df['LOVD_all_clinical'].str.contains('benign')) &
                                 (~filter_df['LOVD_all_clinical'].str.contains('VUS'))) |
                                filter_df['ClinGen_annotation'].str.contains('Pathogenic')]

    # non-actionable heredity variants
    heredity_df = tmp_heredity_df[~tmp_heredity_df.isin(actionable_df.to_dict('l')).all(1)]

    print('------- Heredity variants ----------')
    # print(heredity_df)

    # Uncertain -> LOVD/ClinVar/Clingene not pathogenic(/likely) or not benign(/likely)
    tmp_un_df = filter_df[~filter_df.isin(actionable_df.to_dict('l')).all(1)&
                       ~filter_df.isin(heredity_df.to_dict('l')).all(1)]
    uncertain_df = tmp_un_df[~tmp_un_df['LOVD_all_clinical'].str.contains('benign')&
                             ~tmp_un_df['ClinGen_annotation'].str.contains('Benign')&
                             (~tmp_un_df['CLNSIG'].isin(['Benign', 'Likely_benign']))]

    # print(uncertain_df)
    # uncertain_df.to_csv('uncertain_df.csv', sep=',')

    internal_filter = uncertain_df[(uncertain_df['VAF'].astype(float) >= 0.05) &
                                   (uncertain_df['FAO'].astype(float) >= 10)]
    # ------ COSMIC -------
    COSMIC_df = internal_filter[~internal_filter['cosmic90_coding'].isin(['.'])]
    print('------- COSMIC variants --------')
    # print(COSMIC_df)

    # Prediction
    non_cosmic_df = internal_filter[internal_filter['cosmic90_coding'].isin(['.'])]
    if 'summarized_prediction' in non_cosmic_df.columns:
        unpredict_df = non_cosmic_df[non_cosmic_df['summarized_prediction'] == 'Un_predict']
        tmp_predict_df = non_cosmic_df[~non_cosmic_df.isin(unpredict_df.to_dict('l')).all(1)]
        suspect_df = tmp_predict_df[tmp_predict_df['summarized_prediction'].astype(float) >= 0.71]
        suspect_df = suspect_df.append(unpredict_df)
    else:
        suspect_df = non_cosmic_df[non_cosmic_df['test1$pre_sum'].astype(float) >= 0.71]

    print('------- Suspect variants --------')
    # print(suspect_df)

    return actionable_df, heredity_df, COSMIC_df, suspect_df

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input',
                        help='input source data file name')
    parser.add_argument('--output',
                        help='output path')

    args = parser.parse_args()
    input_file = args.input
    resultFile_path = args.output

    print(input_file.split('\\'))
    tmp = input_file.split('\\')[-1]
    patient_ID = tmp.split('.')[0]
    dirt = '/'.join(input_file.split('\\')[0:-1])
    print(input_file)
    print(dirt)
    print(patient_ID)

    input_df = pd.read_csv(input_file, sep='\t', engine='python')

    actionable_df, heredity_df, COSMIC_df, suspect_df = filter(input_df)
    actionable_num = len(actionable_df)
    heredity_num = len(heredity_df)
    cosmic_num = len(COSMIC_df)
    suspect_num = len(suspect_df)

    act_df_html = actionable_df.to_html(border=0, classes='table table-striped')
    heredity_df_html = heredity_df.to_html(border=0, classes='table table-striped')
    cosmic_df_html = COSMIC_df.to_html(border=0, classes='table table-striped')
    suspect_df_html = suspect_df.to_html(border=0, classes='table table-striped')

    parameters = {'actionable_df': actionable_df, 'actionable_num': actionable_num, 'act_df_html': act_df_html,
                  'heredity_df': heredity_df, 'heredity_num': heredity_num, 'heredity_df_html': heredity_df_html,
                  'COSMIC_df': COSMIC_df, 'cosmic_num': cosmic_num, 'cosmic_df_html': cosmic_df_html,
                  'suspect_df': suspect_df, 'suspect_num': suspect_num, 'suspect_df_html': suspect_df_html,
                  }

    # resultFile_path = output_path + '/' + patient_ID
    #Write pickle file
    print(resultFile_path)

    with open(resultFile_path + '.pickle', 'wb') as wf:
        pickle.dump(parameters, wf)

    #Write Excel file
    writer = pd.ExcelWriter(resultFile_path+ '.xlsx')

    actionable_df.to_excel(writer, sheet_name='Actionable', index=False)
    heredity_df.to_excel(writer, sheet_name='Heredity', index=False)
    COSMIC_df.to_excel(writer, sheet_name='COSMIC', index=False)
    suspect_df.to_excel(writer, sheet_name='Suspect', index=False)

    # Close the Pandas Excel writer and output the Excel file.
    writer.save()

