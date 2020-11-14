#L1000-CTRP-withoutdose-3h dataset
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.metrics import mean_squared_error,r2_score
from sklearn.linear_model import ElasticNet
import warnings
import xgboost as xgb
import scipy as sc
import math
pd.set_option('display.max_columns',15)
warnings.filterwarnings("ignore")
record_start_flag = False
xgboost_model = []
pearson_score = []
mse_score = []

#Retrieving the data
def getData():
    #1.Reading in data
    sig_viability_info = pd.read_table(r"LINCS_CTRP_withoutdose\LINCS_CTRP_withoutdose_3h.txt",sep="	",dtype=str)
    sig_viability_info_filter = sig_viability_info[['signature','cell_line','drug','dose','time','cell_viability']]
    sig_GE_info = pd.read_csv(r"LINCS_CTRP_withoutdose\LINCS_gene_expression_3h.csv",sep=",",dtype=str)
    sig_GE_info_filter = sig_GE_info[:]

    #2.Combining the two datasets
    sig_viability = pd.merge(sig_viability_info_filter,sig_GE_info_filter,left_on='signature',right_on='cid')
    select_sig_viavility2 = pd.DataFrame(sig_viability).loc[:,'780':]
    select_sig_viavility_arr = select_sig_viavility2.copy()
    
    #print(select_sig_viavility_arr)
    #select_sig_viavility1_arr = np.array(select_sig_viavility1)
    #select_sig_viavility2_arr = np.array(select_sig_viavility2)
    #select_sig_viavility_arr = np.insert(select_sig_viavility2_arr,0,values=select_sig_viavility1_arr,axis=1)
    
    #3.Getting the cell viability value
    viability_value = pd.DataFrame(sig_viability).loc[:,'cell_viability']
    
    #4.Dividing the dataset
    xTrain,x_Validation_Test,yTrain,y_Validation_Test = train_test_split(select_sig_viavility_arr,viability_value,test_size=0.3,random_state=50)
    xValidation,xTest,yValidation,yTest = train_test_split(x_Validation_Test,y_Validation_Test,test_size=0.5,random_state=50)
    return xTrain,xValidation,xTest,yTrain,yValidation,yTest

def BinaryToDecimal(x,Pmin,Pmax):
    decimal_number = 0
    for i in range(len(x)):
        decimal_number += x[i]*math.pow(2,i)
    p = Pmin + ((Pmax-Pmin) / (math.pow(2,len(x))-1)) * decimal_number
    return p

#Defining the fitness function
def calc_fitness(binary_code,xSet,ySet,xSet_test,ySet_test):
    #Training the XGBoost model
    dtrain = xgb.DMatrix(xSet,ySet)
    dtest = xgb.DMatrix(xSet_test,ySet_test)
    learning_rate_param = BinaryToDecimal(binary_code[0:6],0.01,0.8)
    n_estimators_param = int(BinaryToDecimal(binary_code[6:12],1000,6000))
    max_depth_param = int(BinaryToDecimal(binary_code[12:18],3,7))
    min_child_weight_param = int(BinaryToDecimal(binary_code[18:24],1,15))
    gamma_param = BinaryToDecimal(binary_code[24:30],0,2)
    subsample_param = BinaryToDecimal(binary_code[30:36],0.01,1)
    colsample_bytree_param = BinaryToDecimal(binary_code[36:42],0,1)
    lambda_param = BinaryToDecimal(binary_code[42:48],0.01,2)
    param = {'silent':True,'obj':'reg:linear','gamma':gamma_param,"eval_metric":"rmse",'max_depth':max_depth_param
         ,'lambda':lambda_param,'subsample':subsample_param,'eta':learning_rate_param,'random_state':50
         ,'colsample_bytree':colsample_bytree_param,'min_child_weight':min_child_weight_param}
    num_round = n_estimators_param
    bst = xgb.train(param,dtrain,num_round)
    resultPredictions = bst.predict(dtest)
    print("Validation_Pearson:",sc.stats.pearsonr(resultPredictions,ySet_test)[0])
    #print("Validation_MSE:",mean_squared_error(resultPredictions,ySet_test))
    if record_start_flag == True:
        xgboost_model.append(bst)
        pearson_score.append(sc.stats.pearsonr(resultPredictions,ySet_test)[0])
        mse_score.append(mean_squared_error(resultPredictions,ySet_test))
    return sc.stats.pearsonr(resultPredictions,ySet_test)[0]

#Obtaining the gene expression data and dividing the datasets into training set and test set
xTrain,xValidation,xTest,yTrain,yValidation,yTest = getData()
#The first step: sorting the importance of the features by using the random forest method and selecting the top features
predictResult_original = []
rfr = RandomForestRegressor(n_estimators=40,max_depth=4)
rfr.fit(xTrain,yTrain)
rfr_pearson = sc.stats.pearsonr(rfr.predict(xValidation),np.array(yValidation,dtype="float32"))[0]
for estimator_index in range(len(rfr.estimators_)):
    predictValue = rfr.estimators_[estimator_index].predict(xValidation)
    predictResult_original.append(sc.stats.pearsonr(predictValue,np.array(yValidation,dtype="float32"))[0])

#Adding random noise to each feature and evaluating their error without noise
np.random.seed(20)
variable_importance = []
for variable_iter in range(xValidation.shape[1]): 
    predictResult_random = []
    random_value = 2 * np.random.rand(xValidation.iloc[:,variable_iter].shape[0]) - 1
    xTest_original = np.array(xValidation.copy(),dtype="float32")
    xTest_original[:,variable_iter] = xTest_original[:,variable_iter] + random_value
    xTest_plus_random = xTest_original
    for estimator_index in range(len(rfr.estimators_)):
        predictValue = rfr.estimators_[estimator_index].predict(xTest_plus_random)
        predictResult_random.append(sc.stats.pearsonr(predictValue,np.array(yValidation,dtype="float32"))[0])
    #Calculating the value of the d_hat
    predict_difference = abs(np.array(predictResult_original) - np.array(predictResult_random))
    d_hat = predict_difference.sum() / len(rfr.estimators_)
    #Calculating the value of the Sd square
    Sd2_sum = 0
    for i in range(len(rfr.estimators_)):
        Sd2_sum += (predict_difference[i] - d_hat)**2
    Sd2 = Sd2_sum / (len(rfr.estimators_) - 1)
    if Sd2 == 0:
        variable_importance.append(0)
    else:
        variable_importance.append(d_hat / math.sqrt(Sd2))
#Sorting the features according to the calculated feature importance
importance_list = {}
importance_list_all = {}
nonzero_flag_index = 0
importance_copy = variable_importance.copy()
importance_sort = np.argsort(importance_copy)[::-1]
for i in range(len(importance_sort)):
    if importance_copy[importance_sort[i]] != 0:
        nonzero_flag_index = i
        importance_list[xValidation.columns[importance_sort[i]]] = i
    else:
        break
#Recording the ranking of all features
for j in range(len(importance_sort)):
    importance_list_all[xValidation.columns[importance_sort[j]]] = j
print('Random Forest Gene Selection OK!')

#The second step: sorting the importance of features by using the elastic network and selecting the top features
validations = []
alphas = []
ratios = []
models = []
alphasrange = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
l1_ratio_range = [0.1,0.2,0.5,0.7,0.7,0.95,1]
#selecting the best parameter by using the validation set
for alpha_index in range(len(alphasrange)):
    for ratio_index in range(len(l1_ratio_range)):
        ElasticNet_Model = ElasticNet(alpha=alphasrange[alpha_index],l1_ratio=l1_ratio_range[ratio_index])
        ElasticNet_Model.fit(xTrain,yTrain)
        models.append(ElasticNet_Model)
        validationResult = ElasticNet_Model.predict(xValidation)
        if math.isnan(sc.stats.pearsonr(validationResult,np.array(yValidation,dtype="float32"))[0]) == False:
            validations.append(sc.stats.pearsonr(validationResult,np.array(yValidation,dtype="float32"))[0])
        else:
            validations.append(0)
        alphas.append(alphasrange[alpha_index])
        ratios.append(l1_ratio_range[ratio_index])
ENet_sortlist_index = np.argsort(validations)[::-1]
ENet_maxpearson_index = ENet_sortlist_index[0]
ENet_pearson = validations[ENet_maxpearson_index]
#print(alphas[ENet_maxpearson_index])
#print(ratios[ENet_maxpearson_index])
#Getting the sequence of characteristic genes
sort_numbers = {}
sort_indexs = np.argsort(abs(models[ENet_maxpearson_index].coef_))[::-1]
for i in range(len(sort_indexs)):
    sort_numbers[xValidation.columns[sort_indexs[i]]] = i
for i in range(len(sort_indexs)):
    if models[0].coef_[sort_indexs[i]] != 0:
        en_nonzero_flag_index = i
    else:
        break

#Calculating the weighted average value
weighted_feature = {}
for weighted_feature_index in range(len(importance_list_all)):
    linear_sort_num = sort_numbers[list(importance_list_all.keys())[weighted_feature_index]]
    weighted_feature[list(importance_list_all.keys())[weighted_feature_index]] = (math.exp(ENet_pearson)*linear_sort_num + math.exp(rfr_pearson)*list(importance_list_all.values())[weighted_feature_index]) / (math.exp(ENet_pearson) + math.exp(rfr_pearson))
#Recording the number of features selected
weighted_feature_nums = (math.exp(ENet_pearson) * en_nonzero_flag_index + math.exp(rfr_pearson) * nonzero_flag_index) / (math.exp(ENet_pearson) + math.exp(rfr_pearson))
selected_genes = sorted(weighted_feature.items(),key=lambda d:d[1],reverse=False)[:int(weighted_feature_nums)]
#print(len(selected_genes))
selected_genes_dict = dict(selected_genes)
print('Elastic Net Gene Selection OK!')

#Dividing the training set and the test set according to the selected features
xTest = np.array(xTest[list(selected_genes_dict.keys())],dtype="float32")
xTrain = np.array(xTrain[list(selected_genes_dict.keys())],dtype="float32")
xValidation = np.array(xValidation[list(selected_genes_dict.keys())],dtype="float32")
yTrain = np.array(yTrain,dtype="float32")
yTest = np.array(yTest,dtype="float32")
yValidation = np.array(yValidation,dtype="float32")

#Using PSO+XGBoost for training
#1.1Particle swarm optimization algorithm initialization
N = 25                 #The number of particles in the group
D = 48                 #The dimension of the particle 
T = 50                  #The maximum number of iterations
c1 = 1.5                 #Learning factor 1
c2 = 1.5                 #Learning factor 2
Wmax = 0.8               #Maximum inertia weight
Wmin = 0.4               #Minimum inertia weight
Vmax = 10                #The maximum speed
Vmin = -10               #The minimum speed
Psai = 2.6

#1.2Initializing the population individuals
#Initializing the binary string
#（Group：learning_rate,n_estimators,max_depth,min_child_weight,gamma,subsample,colsample_bytree）
x = np.random.randint(0,2,(N,D))
v = np.random.rand(N,D) * (Vmax-Vmin) + Vmin
vx = np.random.rand(N,D)

#1.3Initializing the optimal position and optimal value of the individual
record_start_flag = True
record_index = 0
p = x.copy()
pbest = np.ones(N)
for i in range(len(pbest)):
    pbest[i] = calc_fitness(x[i,:],xTrain,yTrain,xValidation,yValidation)

#1.4Initializing the global optimal position and optimal value
g = np.ones(D)
gbest = -float('inf')
for i in range(len(pbest)):
    if pbest[i] > gbest:
        g = p[i,:]
        gbest = pbest[i]
gb = np.ones(T)

#2.Iterating in sequence until the accuracy or number of the iterations has been met
for i in range(T):
    for j in range(N):
        #Updating optimal position and optimal value of the individual
        print('--------------optimize fitness-----------------')
        print('Iterations number:',i,'individual number:',j)
        theFitnessValue = calc_fitness(x[j,:],xTrain,yTrain,xValidation,yValidation)
        record_index = record_index + 1
        if theFitnessValue > pbest[j]:
            p[j,:] = x[j,:]
            pbest[j] = theFitnessValue
        #Updating the global optimal position and optimal value
        #print('pbest[j]>gbest',pbest[j] > gbest)
        if pbest[j] > gbest:
            g = p[j,:]
            gbest = pbest[j]
        #print('g:',g,'T:',i)
        #print('pbest',pbest[j],'gbest',gbest)
        #print('---------------------------')
        #Calculating the dynamic inertia weight value
        IW_alpha_1 = (Wmin*math.exp(Psai) - Wmax*math.exp(2*Psai)) / (1 - math.exp(2*Psai))
        IW_alpha_2 = (Wmax - Wmin*math.exp(Psai)) / (1 - math.exp(2*Psai))
        w = IW_alpha_1 * math.exp(-Psai*i / T) + IW_alpha_2 * math.exp(Psai*i / T)
        #Updating the position and speed values
        v[j,:] = w*v[j,:]+c1*np.random.rand()*(p[j,:]-x[j,:])+c2*np.random.rand()*(g-x[j,:])
        #Boundary condition processing
        for ii in range(D):
            if (v[j,ii] > Vmax) | (v[j,ii] < Vmin):
                v[j,ii] = np.random.rand()*(Vmax-Vmin)+Vmin
        vx[j,:] = 1/(1+np.exp((-1)*np.array(v[j,:])))
        for jj in range(D):
            if vx[j,jj] > np.random.rand():
                x[j,jj] = 1
            else:
                x[j,jj] = 0
    #Recording the global optimal value in the past generations
    gb[i] = gbest
#print(g)
print('gbest:',gbest)

#Evaluating the model by using the test data
best_model_index = pearson_score.index(gbest)
xgboost_best_model = xgboost_model[best_model_index]
dvalidation = xgb.DMatrix(xValidation,yValidation)
dtest = xgb.DMatrix(xTest,yTest)
resultPredictions = xgboost_best_model.predict(dvalidation)
testPredictions = xgboost_best_model.predict(dtest)
print('****************************************')
print('learning_rate_param:',BinaryToDecimal(g[0:6],0.01,0.8))
print('n_estimators_param:',int(BinaryToDecimal(g[6:12],1000,6000)))
print('max_depth_param:',int(BinaryToDecimal(g[12:18],3,7)))
print('min_child_weight_param:',int(BinaryToDecimal(g[18:24],1,15)))
print('gamma_param:',BinaryToDecimal(g[24:30],0,2))
print('subsample_param:',BinaryToDecimal(g[30:36],0.01,1))
print('colsample_bytree_param:',BinaryToDecimal(g[36:42],0,1))
print('lambda_param:',BinaryToDecimal(g[42:48],0.01,2))

print('Validation_Best_Pearson:',sc.stats.pearsonr(resultPredictions,yValidation)[0])
print('Validation_R2:',r2_score(yValidation,resultPredictions))
print("MSE_Validation: %.3f" % (mean_squared_error(resultPredictions,yValidation)))

print('Test_Best_Pearson:',sc.stats.pearsonr(testPredictions,yTest)[0])
print('Test_R2:',r2_score(yTest,testPredictions))
print("MSE_Test: %.3f" % (mean_squared_error(testPredictions,yTest)))

#Saving the model
import pickle
gene_files = open(r"model_adjust\CTRP-L1000-3h_GeneSelection_1","w")
gene_files.write(str(selected_genes_dict))
gene_files.close()
pickle.dump(xgboost_best_model,open(r"model_adjust\CTRP-L1000-3h_XGBoost_1.dat","wb"))