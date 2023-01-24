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
df2 = pd.read_csv(r'./data/Iupred2.csv')
df3 = pd.read_csv(r'./data/PTMinteraction.csv')

df5 = pd.read_csv(r'./data/complex_interaction.csv')
df6 = pd.read_csv(r'./data/complex_interface.csv')

df7 = pd.read_csv(r'./data/complex_atom.csv')
df8=df3[["Organism","Gene","Uniprot","PTM","Site","AA","Int_uniprot","Int_gene","Effect","PMID"]]

df9 = pd.read_csv(r'./data/aa_markpoint.csv')

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


    #基因
    keyword2 = StringField()

    #物种
    search_organism = SelectField('Organism',
                                  choices=[['Human']],
                                  validators=[DataRequired()],
                                  default='Human')

    #PTM
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
    organism = form.search_organism.data
    gene = form.keyword2.data

    AA = form.Phos.data  # Boolean 类型
    BB = form.Ac.data  # Boolean 类型
    CC = form.Me.data  # Boolean 类型
    DD = form.Sumo.data  # Boolean 类型
    EE = form.Ub.data  # Boolean 类型
    FF = form.Glyco.data  # Boolean 类型

    if request.method == 'GET':
        return render_template('Search.html')
    elif request.method == 'POST':
       if str(gene)=="":

           if str(AA) == "True":
               yy1 = df8[(df8.PTM == "Phos")]
           else:
               yy1 = pd.DataFrame( columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect","PMID"])

           if str(BB) == "True":
               yy2 = df8[(df8.PTM == "Ac")]
           else:
               yy2= pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect","PMID"])

           if str(CC) == "True":
               yy3 = df8[(df8.PTM == "Me")]
           else:
               yy3= pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect","PMID"])

           if str(DD) == "True":
               yy4 = df8[(df8.PTM == "Sumo")]
           else:
               yy4= pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect","PMID"])

           if str(EE) == "True":
               yy5 = df8[(df8.PTM == "Ub")]
           else:
               yy5= pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect","PMID"])

           if str(FF) == "True":
               yy6 = df8[(df8.PTM == "Glyco")]
           else:
               yy6= pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect","PMID"])

           yy8 = pd.concat([yy1, yy2, yy3, yy4, yy5, yy6]).drop_duplicates().reset_index(drop=True)

           if len(yy8) == 0:
               yy8= df8[(df8.Organism == organism)].drop_duplicates().reset_index(drop=True)
               yy9 = yy8.to_json(orient='records')
               query_list = json.loads(yy9)
               return render_template('result_test_Anysearch2.html', query=query_list)
           else:
               yy9 = yy8.to_json(orient='records')
               query_list = json.loads(yy9)
               return render_template('result_test_Anysearch2.html', query=query_list)

       else:
           uu1 = df8[df8["Gene"].str.contains(gene, case=False)]
           uu3 = df8[df8["Uniprot"].str.contains(gene, case=False)]
           if len(uu1) == 0:
               uu1 = pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect","PMID"])
           if len(uu3) == 0:
               uu3 = pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect","PMID"])

           uu2=pd.concat([uu1,uu3]).drop_duplicates().reset_index(drop=True)


           if len(uu2) == 0:
               X10 = "Invalid search way!"
               return render_template('Search.html', X10=X10)

           else:
               if str(AA) == "True":
                   yy1 = uu2[(uu2.PTM == "Phos")]
                   if len(yy1)==0:
                       yy1 = pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene","Effect", "PMID"])
               else:
                   yy1 = pd.DataFrame(
                       columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect","PMID"])


               if str(BB) == "True":
                   yy2 = uu2[(uu2.PTM == "Ac")]
                   if len(yy2)==0:
                       yy2 = pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene","Effect", "PMID"])
               else:
                   yy2 = pd.DataFrame(
                       columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect","PMID"])



               if str(CC) == "True":
                   yy3 = uu2[(uu2.PTM == "Me")]
                   if len(yy3)==0:
                       yy3 = pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene","Effect", "PMID"])
               else:
                   yy3 = pd.DataFrame(
                       columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect","PMID"])



               if str(DD) == "True":
                   yy4 = uu2[(uu2.PTM == "Sumo")]
                   if len(yy4)==0:
                       yy4 = pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene","Effect", "PMID"])
               else:
                   yy4 = pd.DataFrame(
                       columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect","PMID"])

               if str(EE) == "True":
                   yy5 = uu2[(uu2.PTM == "Ub")]
                   if len(yy5)==0:
                       yy5 = pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene","Effect", "PMID"])
               else:
                   yy5 = pd.DataFrame(
                       columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect","PMID"])

               if str(FF) == "True":
                   yy6 = uu2[(uu2.PTM == "Glyco")]
                   if len(yy6)==0:
                       yy6 = pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene","Effect", "PMID"])
               else:
                   yy6 = pd.DataFrame(
                       columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect","PMID"])


               yy8 = pd.concat([yy1, yy2, yy3, yy4, yy5, yy6]).drop_duplicates().reset_index(drop=True)

               if len(yy8) == 0:
                  X10 = "Invalid search way!"
                  return render_template('Search.html', X10=X10)
               else:
                  yy9 = yy8.to_json(orient='records')
                  query_list = json.loads(yy9)
                  return render_template('result_test_Anysearch2.html', query=query_list)


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

        # ww5 = df8[df8["Int_uniprot"].str.contains(AAAA, case=False)]
        # if len(ww5) == 0:
        #     ww5 = pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect","PMID"])
        #
        # ww6 = df8[df8["Int_gene"].str.contains(AAAA, case=False)]
        # if len(ww6) == 0:
        #     ww6 = pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect","PMID"])

        ww7 = df8[df8["Effect"].str.contains(AAAA, case=False)]
        if len(ww7) == 0:
            ww7 = pd.DataFrame(columns=["Organism", "Gene", "Uniprot", "PTM", "Site", "AA", "Int_uniprot", "Int_gene", "Effect","PMID"])

        ww8 = pd.concat([ww1, ww2, ww3, ww4,ww7]).drop_duplicates().reset_index(drop=True)

        if len(ww8) == 0:
            X10 = "Invalid search way!"
            return render_template('index.html', X10=X10)
        else:
            ww9 = ww8.to_json(orient='records')
            query_list= json.loads(ww9)
            return render_template('result_test_Anysearch2.html',query=query_list)





# Result page1
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
    
    
    #markpoint
    
    mark=df9[df9.Uniprot ==C1]
    mark1=mark[['Label','Site','Iupred2']]
    mark1.columns=['value','xAxis','yAxis']
    mm2=[{'value': m['value'],'xAxis': m['xAxis'],'yAxis': m['yAxis'],'itemStyle': {'color': "rgba(115, 159, 250, .8)" }} for u,m in mark1.iterrows()]
    


    Y1 = A1 + "_protein_features.csv"
    Y2 = './data/' + Y1
    Iu_hydro.to_csv(Y2, index=False)

    gg = df3[(df3.Uniprot == C1)].reset_index(drop=True)


    # 获得uniprot对应的interaction表格

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
    
    # gg[['Complex','PDBRES']].drop_duplicates()
    
    complexList = gg['Complex'].drop_duplicates().to_list()

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
    K2 = df5[(df5.Complex.isin(complexList))]
    K3 = K2.to_json(orient='records')
    com_int_list = json.loads(K3)
    
    labelList = df7[df7.Complex.isin(complexList)].to_json(orient = "records")


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


    return render_template('result_test_uniprot4.html', A1=A1, B1=B1, C1=C1, D1=D1, E1=E1, F1=F1,
                           I1=I1, J1=J1, K1=K1, Y1=Y1, Z1=Z1, ZZ2=ZZ2, ZZ3=ZZ3,xx=xx,yy=yy,zz=zz,label=label,mm2=mm2,
                           complexList=complexList,
                           labelList=labelList,#[{'Complex':'','x':{{xx}}, 'y':{{yy}}, 'z':{{zz}},'label':"{{label}}",'color':''}]
                           int=int_list,com=com_int_list)

@app.route('/queryFilterdGene/<string:ptmType>')
def queryFilterdGene(ptmType):
    df10=df3[df3.PTM==ptmType]
    gene=df10[['Gene']].drop_duplicates()['Gene'].tolist()
    return jsonify(gene)
    

#Result page2
@app.route('/Search/<string:uniprot>/<string:ptm>/<int:site2>')
def results(uniprot,ptm,site2):
    C1 = uniprot  # uniprot
    proteinptm = ptm
    proteinsite= site2


    basic = df1[(df1.Uniprot == C1)]
    A1 = basic.iat[0, 1]  # 基因名
    B1 = basic.iat[0, 3]  # 蛋白名
    D1 = basic.iat[0, 0]  # 物种名

    Iu_hydro = df2[(df2.Uniprot == C1)]
    E1 = Iu_hydro['Site'].to_list()
    F1 = Iu_hydro['Iupred2'].to_list()

    # markpoint

    mark = df9[(df9.Uniprot == C1)&(df9.PTM==proteinptm)&(df9.Site==proteinsite)]
    print(mark)
    mark1 = mark[['Label', 'Site', 'Iupred2']]
    mark1.columns = ['value', 'xAxis', 'yAxis']
    mm2 = [{'value': m['value'], 'xAxis': m['xAxis'], 'yAxis': m['yAxis'],
            'itemStyle': {'color': "rgba(115, 159, 250, .8)"}} for u, m in mark1.iterrows()]

    Y1 = A1 + "_protein_features.csv"
    Y2 = './data/' + Y1
    Iu_hydro.to_csv(Y2, index=False)

    gg = df3[(df3.Uniprot == C1)& (df3.PTM==proteinptm)&(df3.Site==proteinsite)].reset_index(drop=True)

    # 获得uniprot对应的interaction表格

    gg1 = gg.to_json(orient='records')
    int_list = json.loads(gg1)

    # 获得uniprot对应的int_gene
    ZZ1 = gg[['Int_gene']]
    ZZ1.columns = ['name']
    ZZ1 = ZZ1.drop_duplicates().reset_index(drop=True)

    ZZ2 = [{'id': str(i + 1), 'name': j['name'], 'colors': "#34bf49", "symbolSize": 18} for i, j in ZZ1.iterrows()]
    ZZ3 = [{'source': '0', 'target': str(i + 1)} for i in range(len(ZZ1))]

    # 获得uniprot对应的score
    II = gg[['Site', 'Score']]
    II2 = II.drop_duplicates()
    II3 = II2.sort_values(["Score"], ascending=False)
    I1 = II3['Site'].to_list()
    J1 = II3['Score'].to_list()

    # 获得uniprot对应的复合物信息(first)
    K1 = gg[['Complex']].iat[0, 0]  # 复合物名称

    # gg[['Complex','PDBRES']].drop_duplicates()

    complexList = gg['Complex'].drop_duplicates().to_list()

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
    K2 = df5[(df5.Complex.isin(complexList))]
    K3 = K2.to_json(orient='records')
    com_int_list = json.loads(K3)

    labelList = df7[df7.Complex.isin(complexList)].to_json(orient="records")

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

    return render_template('result_test_uniprot5.html', A1=A1, B1=B1, C1=C1, D1=D1, E1=E1, F1=F1,
                           I1=I1, J1=J1, K1=K1, Y1=Y1, Z1=Z1, ZZ2=ZZ2, ZZ3=ZZ3, xx=xx, yy=yy, zz=zz, label=label,
                           mm2=mm2,
                           complexList=complexList,
                           labelList=labelList,
                           # [{'Complex':'','x':{{xx}}, 'y':{{yy}}, 'z':{{zz}},'label':"{{label}}",'color':''}]
                           int=int_list, com=com_int_list)


#Result page1
@app.route('/Search/<string:ptm>/<string:gene>')
def resultss(ptm,gene):

    proteinptm = ptm
    A1 = gene  # 基因名


    basic = df1[(df1.Gene == A1)]
    
    C1 = basic.iat[0, 2]  #uniprot
    B1 = basic.iat[0, 3]  # 蛋白名
    
    D1 = basic.iat[0, 0]  # 物种名




    Iu_hydro = df2[(df2.Uniprot == C1)]
    E1 = Iu_hydro['Site'].to_list()
    F1 = Iu_hydro['Iupred2'].to_list()

    # markpoint

    mark = df9[(df9.Uniprot == C1)&(df9.PTM==proteinptm)]
    print(mark)
    mark1 = mark[['Label', 'Site', 'Iupred2']]
    mark1.columns = ['value', 'xAxis', 'yAxis']
    mm2 = [{'value': m['value'], 'xAxis': m['xAxis'], 'yAxis': m['yAxis'],
            'itemStyle': {'color': "rgba(115, 159, 250, .8)"}} for u, m in mark1.iterrows()]

    Y1 = A1 + "_protein_features.csv"
    Y2 = './data/' + Y1
    Iu_hydro.to_csv(Y2, index=False)

    gg = df3[(df3.Uniprot == C1)& (df3.PTM==proteinptm)].reset_index(drop=True)

    # 获得uniprot对应的interaction表格

    gg1 = gg.to_json(orient='records')
    int_list = json.loads(gg1)

    # 获得uniprot对应的int_gene
    ZZ1 = gg[['Int_gene']]
    ZZ1.columns = ['name']
    ZZ1 = ZZ1.drop_duplicates().reset_index(drop=True)

    ZZ2 = [{'id': str(i + 1), 'name': j['name'], 'colors': "#34bf49", "symbolSize": 18} for i, j in ZZ1.iterrows()]
    ZZ3 = [{'source': '0', 'target': str(i + 1)} for i in range(len(ZZ1))]

    # 获得uniprot对应的score
    II = gg[['Site', 'Score']]
    II2 = II.drop_duplicates()
    II3 = II2.sort_values(["Score"], ascending=False)
    I1 = II3['Site'].to_list()
    J1 = II3['Score'].to_list()

    # 获得uniprot对应的复合物信息(first)
    K1 = gg[['Complex']].iat[0, 0]  # 复合物名称

    # gg[['Complex','PDBRES']].drop_duplicates()

    complexList = gg['Complex'].drop_duplicates().to_list()

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
    K2 = df5[(df5.Complex.isin(complexList))]
    K3 = K2.to_json(orient='records')
    com_int_list = json.loads(K3)

    labelList = df7[df7.Complex.isin(complexList)].to_json(orient="records")

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

    return render_template('result_test_uniprot6.html', A1=A1, B1=B1, C1=C1, D1=D1, E1=E1, F1=F1,
                           I1=I1, J1=J1, K1=K1, Y1=Y1, Z1=Z1, ZZ2=ZZ2, ZZ3=ZZ3, xx=xx, yy=yy, zz=zz, label=label,
                           mm2=mm2,
                           complexList=complexList,
                           labelList=labelList,
                           # [{'Complex':'','x':{{xx}}, 'y':{{yy}}, 'z':{{zz}},'label':"{{label}}",'color':''}]
                           int=int_list, com=com_int_list)










if __name__ == '__main__':
    # 对应intenal IP
    app.run(host= '127.0.0.1', port=5000, debug=False)

