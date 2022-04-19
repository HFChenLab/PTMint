# -*- coding: utf-8 -*-
"""
Created on March 2022

@author: hill103,Xiaokun Hong
"""

'''
main program for hosting the website PTMint
'''

from flask import Flask, render_template, request, send_from_directory
from flask_bootstrap import Bootstrap
from flask_wtf import FlaskForm
from wtforms import SelectField, SubmitField, RadioField, StringField, BooleanField
from wtforms.validators import DataRequired
import os
import pandas as pd
import zipfile
from flask import jsonify
from flask import jsonify,make_response
import json

COUNT = int(open('static/count.txt').read())
#####################网页相关####################################################
# 载入main csv文件
df1 = pd.read_csv(r'./data/protein_basic_info.csv')
df2 = pd.read_csv(r'./data/Iupred2_Hydropathicity.csv')
df3 = pd.read_csv(r'./data/PTMinteraction.csv')
df4 = pd.read_csv(r'./data/Mutation_info.csv')
df5 = pd.read_csv(r'./data/complex_info.csv')
df6 = pd.read_csv(r'./data/Interface_combine.csv')
df7 = pd.read_csv(r'./data/Complex_atom.csv')
df8=df3[["Organism","Gene","Uniprot","PTM","Site","AA","Int_uniprot","Int_gene","Effect","PMID"]]


class Config:
    SECRET_KEY = 'hard to guess string'
    SSL_DISABLE = False
    WTF_CSRF_ENABLED = False
    DEBUG = False


class NonValidatingSelectField(SelectField):
    # Skip the pre-validation, otherwise it will raise the "Not a valid choice" error
    def pre_validate(self, form):
        pass


class QueryForm(FlaskForm):
    search_type = RadioField('Search Ways',
                             choices=[('uniprot', 'Uniprot'),
                                      ('gene', 'Gene')],
                             validators=[DataRequired()],
                             default='uniprot')
    search_organism = SelectField('Organism',
                                  choices=[['Human', 'Human']],
                                  validators=[DataRequired()],
                                  default='Human')

    uniprot = StringField()
    gene = StringField()

    Phos = BooleanField("Phos", default=True, false_values=('False', 'false', ''))
    Ac = BooleanField("Ac", default=True, false_values=('False', 'false', ''))
    Me = BooleanField("Me", default=True, false_values=('False', 'false', ''))
    Sumo = BooleanField("Sumo", default=True, false_values=('False', 'false', ''))
    Ub = BooleanField("Ub", default=True, false_values=('False', 'false', ''))
    Glyco = BooleanField("Glyco", default=True, false_values=('False', 'false', ''))

    submit = SubmitField('Search')


class QForm(FlaskForm):
    keyword = StringField()
    submit = SubmitField('Search')


#####################主函数####################################################
app = Flask(__name__)
app.config.from_object(Config)
Bootstrap(app)

def global_count():
    global COUNT
    return COUNT
    
app.add_template_global(global_count,'global_count')

# Extra redirect request for "/favicon.ico"
@app.route('/favicon.ico')
def favicon():
    return send_from_directory(os.path.join(app.root_path, 'static'),
                               'favicon.ico', mimetype='image/vnd.microsoft.icon')

@app.before_request
def count():
    global COUNT
    import re
    if not re.match(r'^.*\.\w{1,6}$',request.path):
        COUNT += 1
    if COUNT % 10 ==0:
        open('static/count.txt','w').write(str(COUNT))
    
# Home page
@app.route('/')
def index():
    return render_template('index.html')


# Other pages
@app.route('/<string:page>')
def show(page):
    return render_template('{}.html'.format(page))

# data
@app.route('/data/<string:page>')
def data(page):
    return send_from_directory('data/', page)


# Search database
@app.route('/Search', methods=['GET', 'POST'])
def search():
    form = QueryForm(data=request.form)
    search_type = form.search_type.data
    organism = form.search_organism.data
    AA = form.Phos.data  # Boolean 类型
    BB = form.Ac.data  # Boolean 类型
    CC = form.Me.data  # Boolean 类型
    DD = form.Sumo.data  # Boolean 类型
    EE = form.Ub.data  # Boolean 类型
    FF = form.Glyco.data  # Boolean 类型

    if request.method == 'GET':
        return render_template('Search.html')
    elif request.method == 'POST':
        # 二选一.判断uniprot文本框是否为空.uniprot,(organism+gene)二选一
        if search_type=="uniprot":
            C1 = form.uniprot.data  # uniprot

            if str(AA) == "True":
                g1 = df3[(df3.Uniprot == C1) & (df3.PTM == "Phos")]
            else:
                g1 = pd.DataFrame(
                    columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect",
                             "Method", "Disease", "Co-localization", "PMID", "Complex", "PDBRES", "Gene_chain",
                             "Int_chain", "Score", "Interface", "Domain", "Complex_origin", "PDBID"])

            if str(BB) == "True":
                g2 = df3[(df3.Uniprot == C1) & (df3.PTM == "Ac")]
            else:
                g2 = pd.DataFrame(
                    columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect",
                             "Method", "Disease", "Co-localization", "PMID", "Complex", "PDBRES", "Gene_chain",
                             "Int_chain", "Score", "Interface", "Domain", "Complex_origin", "PDBID"])

            if str(CC) == "True":
                g3 = df3[(df3.Uniprot == C1) & (df3.PTM == "Me")]
            else:
                g3 = pd.DataFrame(
                    columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect",
                             "Method", "Disease", "Co-localization", "PMID", "Complex", "PDBRES", "Gene_chain",
                             "Int_chain", "Score", "Interface", "Domain", "Complex_origin", "PDBID"])

            if str(DD) == "True":
                g4 = df3[(df3.Uniprot == C1) & (df3.PTM == "Sumo")]
            else:
                g4 = pd.DataFrame(
                    columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect",
                             "Method", "Disease", "Co-localization", "PMID", "Complex", "PDBRES", "Gene_chain",
                             "Int_chain", "Score", "Interface", "Domain", "Complex_origin", "PDBID"])

            if str(EE) == "True":
                g5 = df3[(df3.Uniprot == C1) & (df3.PTM == "Ub")]
            else:
                g5 = pd.DataFrame(
                    columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect",
                             "Method", "Disease", "Co-localization", "PMID", "Complex", "PDBRES", "Gene_chain",
                             "Int_chain", "Score", "Interface", "Domain", "Complex_origin", "PDBID"])

            if str(FF) == "True":
                g6 = df3[(df3.Uniprot == C1) & (df3.PTM == "Glyco")]
            else:
                g6 = pd.DataFrame(
                    columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect",
                             "Method", "Disease", "Co-localization", "PMID", "Complex", "PDBRES", "Gene_chain",
                             "Int_chain", "Score", "Interface", "Domain", "Complex_origin", "PDBID"])

            gg = pd.concat([g1, g2, g3, g4, g5, g6]).reset_index(drop=True)

            basic = df1[(df1.Uniprot == C1)]

            if len(gg) == 0:
                X1 = "Invalid search way!"
                return render_template('Search.html', X1=X1)
            else:
                A1 = basic.iat[0, 1]  # 基因名
                B1 = basic.iat[0, 3]  # 蛋白名
                D1 = basic.iat[0, 0]  # 物种名

                Iu_hydro = df2[(df2.Uniprot == C1)]
                E1 = Iu_hydro['Site'].to_list()
                F1 = Iu_hydro['Iupred2'].to_list()
                G1 = Iu_hydro['Site'].to_list()
                H1 = Iu_hydro['Hydropathicity'].to_list()

                Y1= A1 + "_protein_features.csv"
                Y2= './data/' +Y1
                Iu_hydro.to_csv(Y2, index=False)




                # 获得uniprot对应的interaction表格
                # gg1 = gg.to_json(orient="records")
                # gg2 = gg1.replace("[", '{\n"data":[')
                # gg3 = gg2.replace("]", "]\n}")
                # with open("./data/interaction.txt", "w") as f:
                #     f.write(gg3)

                gg1 = gg.to_json(orient='records')
                int_list = json.loads(gg1)

                # 获得uniprot对应的mutation表格
                if str(AA) == "True":
                    h1 = df4[(df4.Uniprot == C1) & (df4.PTM == "Phos")]
                else:
                    h1 = pd.DataFrame(
                        columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "GenBank_ID", "Mutation_position",
                                 "Mutation_alt", "Mutation_type"])

                if str(BB) == "True":
                    h2 = df4[(df4.Uniprot == C1) & (df4.PTM == "Ac")]
                else:
                    h2 = pd.DataFrame(
                        columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "GenBank_ID", "Mutation_position",
                                 "Mutation_alt", "Mutation_type"])

                if str(CC) == "True":
                    h3 = df4[(df4.Uniprot == C1) & (df4.PTM == "Me")]
                else:
                    h3 = pd.DataFrame(
                        columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "GenBank_ID", "Mutation_position",
                                 "Mutation_alt", "Mutation_type"])

                if str(DD) == "True":
                    h4 = df4[(df4.Uniprot == C1) & (df4.PTM == "Sumo")]
                else:
                    h4 = pd.DataFrame(
                        columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "GenBank_ID", "Mutation_position",
                                 "Mutation_alt", "Mutation_type"])

                if str(EE) == "True":
                    h5 = df4[(df4.Uniprot == C1) & (df4.PTM == "Ub")]
                else:
                    h5 = pd.DataFrame(
                        columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "GenBank_ID", "Mutation_position",
                                 "Mutation_alt", "Mutation_type"])

                if str(FF) == "True":
                    h6 = df4[(df4.Uniprot == C1) & (df4.PTM == "Glyco")]
                else:
                    h6 = pd.DataFrame(
                        columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "GenBank_ID", "Mutation_position",
                                 "Mutation_alt", "Mutation_type"])

                hh = pd.concat([h1, h2, h3, h4, h5, h6]).reset_index(drop=True)

                hh1 = hh.to_json(orient="records")
                hh2 = hh1.replace("[", '{\n"data":[')
                hh3 = hh2.replace("]", "]\n}")
                with open("./data/mutation.txt", "w") as f:
                    f.write(hh3)

                #获得uniprot对应的int_gene
                ZZ1 = gg[['Int_gene']]
                ZZ1.columns = ['name']
                ZZ1 = ZZ1.drop_duplicates().reset_index(drop=True)

                ZZ2 = [{'id': str(i + 1), 'name': j['name'],'colors':"#34bf49","symbolSize": 18} for i, j in ZZ1.iterrows()]
                ZZ3 = [{'source': '0', 'target': str(i+1)} for i in range(len(ZZ1))]



                # 获得uniprot对应的score
                II = gg[['Site', 'Score']]
                II2 = II.drop_duplicates()
                II3 = II2.sort_values(["Score"], ascending=False)
                I1 = II3['Site'].to_list()
                J1 = II3['Score'].to_list()

                # 获得uniprot对应的复合物信息(first)
                K1 = gg[['Complex']].iat[0, 0]  # 复合物名称

                # 获得uniprot对应的复合物信息(first)对应的pdbres
                U1 = gg[['PDBRES']].iat[0, 0]  # pdbres
                U2 = gg[['PTM']].iat[0, 0]  # PTM
                U3 = gg[['Site']].iat[0, 0]  # Site

                U4 = df7[(df7.Complex == K1) & (df7.PDBRES == U1) & (df7.PTM == U2) & (df7.Site == U3)]
                xx = U4.iat[0, 6]
                yy = U4.iat[0, 7]
                zz = U4.iat[0, 8]
                label = U4.iat[0, 4]



                #获得第一个复合物interfaction的信息
                # K2 = df5[(df5.Complex == K1)]
                # K3 = K2.to_json(orient="records")
                # K4 = K3.replace("[", '{\n"data":[')
                # K5 = K4.replace("]", "]\n}")
                # with open("./data/complex.txt", "w") as f:
                #     f.write(K5)

                K2 = df5[(df5.Complex == K1)]
                K3 = K2.to_json(orient='records')
                com_int_list = json.loads(K3)

                # 获取uniprot对应的复合物(all)并压缩成一个文件
                # 定义压缩的函数
                def compress_attaches(List, out_name):
                    f = zipfile.ZipFile(out_name, 'w', zipfile.ZIP_DEFLATED)
                    for file in List:
                        f.write(file)
                    f.close()

                P = gg[['Complex']].drop_duplicates()
                P['path'] = "./data/"
                P['pdb'] = ".pdb"
                P['concat'] = P['path'] + P['Complex'] + P['pdb']
                List1 = P['concat'].tolist()
                name = P['Complex'].tolist()

                #得到com_int.csv文件
                interact = pd.DataFrame(columns=['Complex', 'AA1', 'AA2', 'Type'])
                for e in name:
                    interact1 = df5[(df5.Complex == e)]
                    interact = pd.concat([interact,interact1])
                interact.to_csv('./data/com_int.csv',index=False)

                # 得到com_interface.csv文件
                interface = pd.DataFrame(columns=['Complex', 'Chain', 'Site'])
                for f in name:
                    interface1 = df6[(df6.Complex == f)]
                    interface = pd.concat([interface,interface1])
                interface.to_csv('./data/com_interface.csv',index=False)

                List2 = ['./data/com_int.csv', './data/com_interface.csv']
                List = List1 + List2

                #将pdb文件和两个csv文件压缩
                Z1 = A1 + "_complex_analysis.zip"
                Z2='./data/'+ Z1

                compress_attaches(List, Z2)

                return render_template('result_test_uniprot.html', A1=A1, B1=B1, C1=C1, D1=D1, E1=E1, F1=F1, G1=G1, H1=H1,
                                       I1=I1, J1=J1, K1=K1,Y1=Y1,Z1=Z1,xx=xx,yy=yy,zz=zz,ZZ2=ZZ2,ZZ3=ZZ3,
                                       int=int_list,com=com_int_list,label=label)

        else:

            A1 = form.gene.data  # 基因名
            D1 = organism  # 物种名

            if str(AA) == "True":
                g1 = df3[(df3.Gene == A1) & (df3.Organism == D1) & (df3.PTM == "Phos")]
            else:
                g1 = pd.DataFrame(
                    columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect",
                             "Method", "Disease", "Co-localization", "PMID", "Complex", "PDBRES", "Gene_chain",
                             "Int_chain", "Score", "Interface", "Domain", "Complex_origin", "PDBID"])

            if str(BB) == "True":
                g2 = df3[(df3.Gene == A1) & (df3.Organism == D1) & (df3.PTM == "Ac")]
            else:
                g2 = pd.DataFrame(
                    columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect",
                             "Method", "Disease", "Co-localization", "PMID", "Complex", "PDBRES", "Gene_chain",
                             "Int_chain", "Score", "Interface", "Domain", "Complex_origin", "PDBID"])

            if str(CC) == "True":
                g3 = df3[(df3.Gene == A1) & (df3.Organism == D1) & (df3.PTM == "Me")]
            else:
                g3 = pd.DataFrame(
                    columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect",
                             "Method", "Disease", "Co-localization", "PMID", "Complex", "PDBRES", "Gene_chain",
                             "Int_chain", "Score", "Interface", "Domain", "Complex_origin", "PDBID"])

            if str(DD) == "True":
                g4 = df3[(df3.Gene == A1) & (df3.Organism == D1) & (df3.PTM == "Sumo")]
            else:
                g4 = pd.DataFrame(
                    columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect",
                             "Method", "Disease", "Co-localization", "PMID", "Complex", "PDBRES", "Gene_chain",
                             "Int_chain", "Score", "Interface", "Domain", "Complex_origin", "PDBID"])

            if str(EE) == "True":
                g5 = df3[(df3.Gene == A1) & (df3.Organism == D1) & (df3.PTM == "Ub")]
            else:
                g5 = pd.DataFrame(
                    columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect",
                             "Method", "Disease", "Co-localization", "PMID", "Complex", "PDBRES", "Gene_chain",
                             "Int_chain", "Score", "Interface", "Domain", "Complex_origin", "PDBID"])

            if str(FF) == "True":
                g6 = df3[(df3.Gene == A1) & (df3.Organism == D1) & (df3.PTM == "Glyco")]
            else:
                g6 = pd.DataFrame(
                    columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect",
                             "Method", "Disease", "Co-localization", "PMID", "Complex", "PDBRES", "Gene_chain",
                             "Int_chain", "Score", "Interface", "Domain", "Complex_origin", "PDBID"])

            gg = pd.concat([g1, g2, g3, g4, g5, g6]).reset_index(drop=True)

            basic = df1[(df1.Gene == A1) & (df1.Organism == D1)]

            if len(gg) == 0:
                X1 = "Invalid search way!"
                return render_template('Search.html', X1=X1)
            else:
                B1 = basic.iat[0, 3]  # 蛋白名
                C1 = basic.iat[0, 2]  # Uniprot
                Iu_hydro = df2[(df2.Uniprot == C1)]
                E1 = Iu_hydro['Site'].to_list()
                F1 = Iu_hydro['Iupred2'].to_list()
                F1 = Iu_hydro['Iupred2'].to_list()
                G1 = Iu_hydro['Site'].to_list()
                H1 = Iu_hydro['Hydropathicity'].to_list()

                Y1= A1 + "_protein_features.csv"
                Y2= './data/' +Y1
                Iu_hydro.to_csv(Y2, index=False)

                # 获得uniprot对应的interaction表格
                # gg1 = gg.to_json(orient="records")
                # gg2 = gg1.replace("[", '{\n"data":[')
                # gg3 = gg2.replace("]", "]\n}")
                # with open("./data/interaction.txt", "w") as f:
                #     f.write(gg3)
                gg1 = gg.to_json(orient='records')
                int_list = json.loads(gg1)

                # 获得uniprot对应的mutation表格
                if str(AA) == "True":
                    h1 = df4[(df4.Uniprot == C1) & (df4.PTM == "Phos")]
                else:
                    h1 = pd.DataFrame(
                        columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "GenBank_ID", "Mutation_position",
                                 "Mutation_alt", "Mutation_type"])

                if str(BB) == "True":
                    h2 = df4[(df4.Uniprot == C1) & (df4.PTM == "Ac")]
                else:
                    h2 = pd.DataFrame(
                        columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "GenBank_ID", "Mutation_position",
                                 "Mutation_alt", "Mutation_type"])

                if str(CC) == "True":
                    h3 = df4[(df4.Uniprot == C1) & (df4.PTM == "Me")]
                else:
                    h3 = pd.DataFrame(
                        columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "GenBank_ID", "Mutation_position",
                                 "Mutation_alt", "Mutation_type"])

                if str(DD) == "True":
                    h4 = df4[(df4.Uniprot == C1) & (df4.PTM == "Sumo")]
                else:
                    h4 = pd.DataFrame(
                        columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "GenBank_ID", "Mutation_position",
                                 "Mutation_alt", "Mutation_type"])

                if str(EE) == "True":
                    h5 = df4[(df4.Uniprot == C1) & (df4.PTM == "Ub")]
                else:
                    h5 = pd.DataFrame(
                        columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "GenBank_ID", "Mutation_position",
                                 "Mutation_alt", "Mutation_type"])

                if str(FF) == "True":
                    h6 = df4[(df4.Uniprot == C1) & (df4.PTM == "Glyco")]
                else:
                    h6 = pd.DataFrame(
                        columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "GenBank_ID", "Mutation_position",
                                 "Mutation_alt", "Mutation_type"])

                hh = pd.concat([h1, h2, h3, h4, h5, h6]).reset_index(drop=True)

                hh1 = hh.to_json(orient="records")
                hh2 = hh1.replace("[", '{\n"data":[')
                hh3 = hh2.replace("]", "]\n}")
                with open("./data/mutation.txt", "w") as f:
                    f.write(hh3)

                #获得uniprot对应的int_gene
                ZZ1 = gg[['Int_gene']]
                ZZ1.columns = ['name']
                ZZ1 = ZZ1.drop_duplicates().reset_index(drop=True)

                ZZ2 = [{'id': str(i + 1), 'name': j['name'],'colors':"#34bf49","symbolSize": 18} for i, j in ZZ1.iterrows()]
                ZZ3 = [{'source': '0', 'target': str(i+1)} for i in range(len(ZZ1))]

                # 获得uniprot对应的score
                II = gg[['Site', 'Score']]
                II2 = II.drop_duplicates()
                II3 = II2.sort_values(["Score"], ascending=False)
                I1 = II3['Site'].to_list()
                J1 = II3['Score'].to_list()

                # 获得uniprot对应的复合物信息first
                K1 = gg[['Complex']].iat[0, 0]  # 复合物名称

                # 获得uniprot对应的复合物信息(first)对应的pdbres
                U1 = gg[['PDBRES']].iat[0, 0]  # pdbres
                U2 = gg[['PTM']].iat[0, 0]  # PTM
                U3 = gg[['Site']].iat[0, 0]  # Site

                U4 = df7[(df7.Complex == K1) & (df7.PDBRES == U1) & (df7.PTM == U2) & (df7.Site == U3)]
                xx = U4.iat[0, 6]
                yy = U4.iat[0, 7]
                zz = U4.iat[0, 8]
                label = U4.iat[0, 4]

                # 获得第一个复合物interfaction的信息
                # K2 = df5[(df5.Complex == K1)]
                # K3 = K2.to_json(orient="records")
                # K4 = K3.replace("[", '{\n"data":[')
                # K5 = K4.replace("]", "]\n}")
                # with open("./data/complex.txt", "w") as f:
                #     f.write(K5)

                K2 = df5[(df5.Complex == K1)]
                K3 = K2.to_json(orient='records')
                com_int_list = json.loads(K3)

                # 获取uniprot对应的复合物(all)并压缩成一个文件
                # 定义压缩的函数
                def compress_attaches(List, out_name):
                    f = zipfile.ZipFile(out_name, 'w', zipfile.ZIP_DEFLATED)
                    for file in List:
                        f.write(file)
                    f.close()

                P = gg[['Complex']].drop_duplicates()
                P['path'] = "./data/"
                P['pdb'] = ".pdb"
                P['concat'] = P['path'] + P['Complex'] + P['pdb']
                List1 = P['concat'].tolist()
                name = P['Complex'].tolist()

                # 得到com_int.csv文件
                interact = pd.DataFrame(columns=['Complex', 'AA1', 'AA2', 'Type'])
                for e in name:
                    interact1 = df5[(df5.Complex == e)]
                    interact = interact.append(interact1)
                interact.to_csv('./data/com_int.csv',index=False)

                # 得到com_interface.csv文件
                interface = pd.DataFrame(columns=['Complex', 'Chain', 'Site'])
                for f in name:
                    interface1 = df6[(df6.Complex == f)]
                    interface = interface.append(interface1)
                interface.to_csv('./data/com_interface.csv',index=False)

                List2 = ['./data/com_int.csv', './data/com_interface.csv']
                List = List1 + List2

                #将pdb文件和两个csv文件压缩
                Z1 = A1 + "_complex_analysis.zip"
                Z2='./data/'+ Z1

                compress_attaches(List, Z2)

                return render_template('result_test_uniprot.html', A1=A1, B1=B1, C1=C1, D1=D1, E1=E1, F1=F1, G1=G1, H1=H1,
                                       I1=I1, J1=J1, K1=K1,Y1=Y1,Z1=Z1,xx=xx,yy=yy,zz=zz,ZZ2=ZZ2,ZZ3=ZZ3,
                                       int=int_list,com=com_int_list,label=label)


@app.route('/Citation', methods=['GET', 'POST'])
def citation():
    form = QForm(data=request.form)
    AAAA = form.keyword.data

    if request.method == 'GET':
        return render_template('Citation.html')
    elif request.method == 'POST':
        ww1 = df1[df1["Organism"].str.contains(AAAA, case=False)]
        if len(ww1) == 0:
            ww1 = pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "Protein_name", "Protein_length"])

        ww2 = df1[df1["Gene"].str.contains(AAAA, case=False)]
        if len(ww2) == 0:
            ww2 = pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "Protein_name", "Protein_length"])

        ww3 = df1[df1["Uniprot"].str.contains(AAAA, case=False)]
        if len(ww3) == 0:
            ww3 = pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "Protein_name", "Protein_length"])

        ww4 = df1[df1["Protein_name"].str.contains(AAAA, case=False)]
        if len(ww4) == 0:
            ww4 = pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "Protein_name", "Protein_length"])

        ww5 = pd.concat([ww1, ww2, ww3, ww4]).drop_duplicates().reset_index(drop=True)

        if len(ww5) == 0:
            X10 = "Invalid search way!"
            return render_template('Citation.html', X10=X10)
        else:
            ww6 = ww5.to_json(orient="records")
            ww7 = ww6.replace("[", '{\n"data":[')
            ww8 = ww7.replace("]", "]\n}")
            with open("./data/search.txt", "w") as f:
                f.write(ww8)
            return render_template('result_test_Anysearch.html')

# @app.route('/Home', methods=['GET', 'POST'])
# def index2():
#     form = QForm(data=request.form)
#     AAAA = form.keyword.data
#
#     if request.method == 'GET':
#         return render_template('index.html')
#     elif request.method == 'POST':
#         ww1 = df1[df1["Organism"].str.contains(AAAA, case=False)]
#         if len(ww1) == 0:
#             ww1 = pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "Protein_name", "Protein_length"])
#
#         ww2 = df1[df1["Gene"].str.contains(AAAA, case=False)]
#         if len(ww2) == 0:
#             ww2 = pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "Protein_name", "Protein_length"])
#
#         ww3 = df1[df1["Uniprot"].str.contains(AAAA, case=False)]
#         if len(ww3) == 0:
#             ww3 = pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "Protein_name", "Protein_length"])
#
#         ww4 = df1[df1["Protein_name"].str.contains(AAAA, case=False)]
#         if len(ww4) == 0:
#             ww4 = pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "Protein_name", "Protein_length"])
#
#         ww5 = pd.concat([ww1, ww2, ww3, ww4]).drop_duplicates().reset_index(drop=True)
#
#         if len(ww5) == 0:
#             X10 = "Invalid search way!"
#             return render_template('index.html', X10=X10)
#         else:
#             ww6 = ww5.to_json(orient="records")
#             ww7 = ww6.replace("[", '{\n"data":[')
#             ww8 = ww7.replace("]", "]\n}")
#             with open("./data/search.txt", "w") as f:
#                 f.write(ww8)
#             return render_template('result_test_Anysearch.html')


@app.route('/Home', methods=['GET', 'POST'])
def index2():
    form = QForm(data=request.form)
    AAAA = form.keyword.data

    if request.method == 'GET':
        return render_template('index.html')
    elif request.method == 'POST':
        ww1 = df8[df8["Organism"].str.contains(AAAA, case=False)]
        if len(ww1) == 0:
            ww1 = pd.DataFrame(columns=["Organism","Gene","Uniprot","PTM","Site","AA","Int_uniprot","Int_gene","Effect","PMID"])

        ww2 = df8[df8["Gene"].str.contains(AAAA, case=False)]
        if len(ww2) == 0:
            ww2 = pd.DataFrame(columns=["Organism","Gene","Uniprot","PTM","Site","AA","Int_uniprot","Int_gene","Effect","PMID"])

        ww3 = df8[df8["Uniprot"].str.contains(AAAA, case=False)]
        if len(ww3) == 0:
            ww3 = pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect","PMID"])

        ww4 = df8[df8["PTM"].str.contains(AAAA, case=False)]
        if len(ww4) == 0:
            ww4 = pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect","PMID"])

        ww5 = df8[df8["Int_uniprot"].str.contains(AAAA, case=False)]
        if len(ww5) == 0:
            ww5 = pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect","PMID"])

        ww6 = df8[df8["Int_gene"].str.contains(AAAA, case=False)]
        if len(ww6) == 0:
            ww6 = pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect","PMID"])

        ww7 = df8[df8["Effect"].str.contains(AAAA, case=False)]
        if len(ww7) == 0:
            ww7 = pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect","PMID"])

        ww8 = pd.concat([ww1, ww2, ww3, ww4,ww5,ww6,ww7]).drop_duplicates().reset_index(drop=True)

        if len(ww8) == 0:
            X10 = "Invalid search way!"
            return render_template('index.html', X10=X10)
        else:
            ww9 = ww8.to_json(orient='records')
            query_list= json.loads(ww9)
            print(type(query_list))
            return render_template('result_test_Anysearch2.html',query=query_list)



# Result page
@app.route('/Search/<string:uniprot>')
def result(uniprot):
    C1=uniprot  #uniprot
    basic = df1[(df1.Uniprot == C1)]
    A1 = basic.iat[0, 1]  # 基因名
    B1 = basic.iat[0, 3]  # 蛋白名
    D1 = basic.iat[0, 0]  # 物种名

    Iu_hydro = df2[(df2.Uniprot == C1)]
    E1 = Iu_hydro['Site'].to_list()
    F1 = Iu_hydro['Iupred2'].to_list()
    G1 = Iu_hydro['Site'].to_list()
    H1 = Iu_hydro['Hydropathicity'].to_list()

    Y1 = A1 + "_protein_features.csv"
    Y2 = './data/' + Y1
    Iu_hydro.to_csv(Y2, index=False)

    gg = df3[(df3.Uniprot == C1)].reset_index(drop=True)


    # 获得uniprot对应的interaction表格
    # gg1 = gg.to_json(orient="records")
    # gg2 = gg1.replace("[", '{\n"data":[')
    # gg3 = gg2.replace("]", "]\n}")
    # with open("./data/interaction.txt", "w") as f:
    #     f.write(gg3)
    gg1 = gg.to_json(orient='records')
    int_list = json.loads(gg1)



    # 获得uniprot对应的int_gene
    ZZ1 = gg[['Int_gene']]
    ZZ1.columns = ['name']
    ZZ1 = ZZ1.drop_duplicates().reset_index(drop=True)

    ZZ2 = [{'id': str(i + 1), 'name': j['name'],'colors':"#34bf49","symbolSize": 18} for i, j in ZZ1.iterrows()]
    ZZ3 = [{'source': '0', 'target': str(i + 1)} for i in range(len(ZZ1))]

    # 获得uniprot对应的score
    II = gg[['Site', 'Score']]
    II2 = II.drop_duplicates()
    II3 = II2.sort_values(["Score"], ascending=False)
    I1 = II3['Site'].to_list()
    J1 = II3['Score'].to_list()

    # 获得uniprot对应的复合物信息(first)
    K1 = gg[['Complex']].iat[0, 0]  # 复合物名称

    # 获得uniprot对应的复合物信息(first)对应的pdbres
    U1 = gg[['PDBRES']].iat[0, 0]  # pdbres
    U2=  gg[['PTM']].iat[0, 0]  #PTM
    U3 = gg[['Site']].iat[0, 0]  # Site

    U4 = df7[(df7.Complex == K1) & (df7.PDBRES == U1) & (df7.PTM == U2) & (df7.Site == U3)]
    xx = U4.iat[0, 6]
    yy = U4.iat[0, 7]
    zz = U4.iat[0, 8]
    label=U4.iat[0, 4]

    # 获得第一个复合物interfaction的信息
    # K2 = df5[(df5.Complex == K1)]
    # K3 = K2.to_json(orient="records")
    # K4 = K3.replace("[", '{\n"data":[')
    # K5 = K4.replace("]", "]\n}")
    # with open("./data/complex.txt", "w") as f:
    #     f.write(K5)

    K2 = df5[(df5.Complex == K1)]
    K3 = K2.to_json(orient='records')
    com_int_list = json.loads(K3)


    # 获取uniprot对应的复合物(all)并压缩成一个文件
    # 定义压缩的函数
    def compress_attaches(List, out_name):
        f = zipfile.ZipFile(out_name, 'w', zipfile.ZIP_DEFLATED)
        for file in List:
            f.write(file)
        f.close()

    P = gg[['Complex']].drop_duplicates()
    P['path'] = "./data/"
    P['pdb'] = ".pdb"
    P['concat'] = P['path'] + P['Complex'] + P['pdb']
    List1 = P['concat'].tolist()
    name = P['Complex'].tolist()

    # 得到com_int.csv文件
    interact = pd.DataFrame(columns=['Complex', 'AA1', 'AA2', 'Type'])
    for e in name:
        interact1 = df5[(df5.Complex == e)]
        interact = pd.concat([interact, interact1])
    interact.to_csv('./data/com_int.csv', index=False)

    # 得到com_interface.csv文件
    interface = pd.DataFrame(columns=['Complex', 'Chain', 'Site'])
    for f in name:
        interface1 = df6[(df6.Complex == f)]
        interface = pd.concat([interface, interface1])
    interface.to_csv('./data/com_interface.csv', index=False)

    List2 = ['./data/com_int.csv', './data/com_interface.csv']
    List = List1 + List2

    # 将pdb文件和两个csv文件压缩
    Z1 = A1 + "_complex_analysis.zip"
    Z2 = './data/' + Z1

    compress_attaches(List, Z2)


    return render_template('result_test_uniprot2.html', A1=A1, B1=B1, C1=C1, D1=D1, E1=E1, F1=F1, G1=G1, H1=H1,
                           I1=I1, J1=J1, K1=K1, Y1=Y1, Z1=Z1, xx=xx, yy=yy, zz=zz, ZZ2=ZZ2, ZZ3=ZZ3,
                           int=int_list,com=com_int_list,label=label)



if __name__ == '__main__':
    # 对应intenal IP
    app.run(host='192.168.1.201', port=8082, debug=False)


