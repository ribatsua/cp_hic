from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.model_selection import train_test_split
import numpy as np
import pandas as pd
import shap
import h5py
from tqdm import tqdm
from sklearn.metrics import classification_report
import sys

from data_label_process import process_data, load_dataset, load_data_to_shap
from cnn_utils import  Simple1DCNN, train,multi_input_cnn
from keras.utils import plot_model  
 

def load_sim_data(data,label,chr_index):
    data_np=np.array(data)
    label_encoder = LabelEncoder()
    indexed_labels = label_encoder.fit_transform(label)
    onehot_encoder = OneHotEncoder(sparse=False)
    encoded_labels = onehot_encoder.fit_transform(indexed_labels.reshape(-1, 1))
    encoded_labels_df = pd.DataFrame(encoded_labels, columns=[f'_{i}' for i in range(encoded_labels.shape[1])])
    X_train, X_test, y_train, y_test = train_test_split(data_np.transpose(), encoded_labels_df, test_size=0.2, random_state=123)
    X_train_3d = np.expand_dims(X_train, axis=2)
    X_test_3d = np.expand_dims(X_test, axis=2)
    print("X_train.shape:", X_train_3d.shape)
    print("y_train.shape:", y_train.shape)
    print()
    print("X_test.shape:", X_test_3d.shape)
    print("y_test.shape:", y_test.shape)
    # train_dataset =[]
    # test_dataset = []
    # input_shapes_list =[]
    # for i in range(len(chr_index)):
    # # for i in range(2):
    #     train_dataset.append(X_train_3d[:,chr_index[i],:])
    #     input_shapes_list.append(X_train_3d[:,chr_index[i],:].shape[1:])
    #     test_dataset.append(X_test_3d[:,chr_index[i],:])

    dataset = (X_train_3d,y_train,X_test_3d,y_test) 

    return dataset

def run_shap(model, data ,label,parm,chr_index):
    group='compartment_raw'
    data_np=np.array(data)
    label_encoder = LabelEncoder()
    indexed_labels = label_encoder.fit_transform(label)
    onehot_encoder = OneHotEncoder(sparse=False)
    encoded_labels = onehot_encoder.fit_transform(indexed_labels.reshape(-1, 1))

    encoded_labels_df = pd.DataFrame(encoded_labels, columns=[f'_{i}' for i in range(encoded_labels.shape[1])])
    X = np.expand_dims(data_np.transpose(), axis=2)
    # temp_data = np.expand_dims(data_np.transpose(), axis=2)
    # print(temp_data.shape)
    # X =[]
   
    # for i in range(len(chr_index)):

    #     X.append(temp_data[:,chr_index[i],:])
    y = np.expand_dims(encoded_labels_df,axis=2)
    predictions = model.predict(X, batch_size=16)
    print(classification_report(y.argmax(axis=1),
    predictions.argmax(axis=1), target_names=label_encoder.classes_))

    # background = [array[np.random.choice(array.shape[0], 100, replace=False)] for array in X]
    # explainer = shap.GradientExplainer(model, background)
    # out_list = []
    # num_samples = X[1].shape[0]
   
    # for sample in tqdm(range(num_samples)):
       
    #     shap_values = explainer.shap_values([array[sample : sample + 1] for array in X])
    #     out_list.append(shap_values)
    # concatenated_third_level = [[np.concatenate(array_list, axis=1) for array_list in sublist] for sublist in out_list]
    # shap_arr = np.squeeze(np.array(concatenated_third_level))
    background=X[np.random.choice(X.shape[0], 100, replace=False)]
    e = shap.DeepExplainer(model, background)
    print("e.expected_value:",e.expected_value)
    out_list = []
    num_samples = np.shape(X)[0]
    for sample in tqdm(range(num_samples)):
        shap_values = e.shap_values(X[sample : sample + 1])
        out_list.append(shap_values)
    shap_arr = np.squeeze(np.array(out_list)) 


    # file_path='/home/dataset/sim_data/cellcycle/1m/SHAP/'+parm+'_shap.h5'
    file_path='/home/dataset/sim_data/cellcycle/250k/SHAP/'+parm+'_shap.h5'
    with h5py.File(file_path, 'w') as hf:
        grp = hf.create_group(group)
        for i in range(shap_arr.shape[2]):
            # cell_data = shap_arr[i].T  
            bin_data=shap_arr[:,:,i] 
            # grp.create_dataset(f'cell_{i + 1}', data=cell_data)
            grp.create_dataset(f'bin_{i + 1}', data=bin_data)
        label_shap=grp.create_group('label')
        label_shap.create_dataset('cell_type', data=[l.encode('utf8') for l in label],
                            dtype=h5py.special_dtype(vlen=str))
    
    hf.close()
    print("SHAP_value save to:",file_path)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("用法: python script.py <参数>")
        sys.exit(1)

    parm = sys.argv[1]
    # 使用参数进行处理

    # parm = 'var4'
    print(f"处理参数: {parm}")
    # data1_path = '/home/sim_data/meanA-1.csv'
    # data2_path = '/home/sim_data/var1.5.csv'

    # data1_path = '/home/dataset/sim_data/cellcycle/1m/compartment_raw_0_1171.csv_mean1+var1.csv'
    # data2_path = '/home/dataset/sim_data/cellcycle/1m/compartment_raw_0_1171.csv_'+parm+'.csv'

    chr_list = ["chr1","chr2","chr3","chr4","chr5",
    "chr6","chr7","chr8","chr9","chr10",
    "chr11","chr12","chr13","chr14","chr15",
    "chr16","chr17","chr18","chr19",'chrX']
    # X=h5py.File('/home/dataset/cellcycle/1m/tmp/scCompartment_cellcycle.hdf5', 'r')
    X=h5py.File('/home/dataset/cellcycle/sparse_input/250k/raw/compartment250k.hdf5', 'r')
    sim_data=h5py.File('/home/dataset/sim_data/cellcycle/250k/compartment250k_simulate.hdf5','r')

    bin_chr = np.array([chrom.decode('utf-8') for chrom in X['compartment_raw']['bin']['chrom']])
    chr_index ={}
        # chr_data =[]
    print(chr_list)
    for i in range(len(chr_list)):
            chr_index[i] = np.where(bin_chr == chr_list[i])[0]
    
    cell_keys = list(filter(lambda key: 'cell_' in key, sim_data['compartment250k_mean1+var1'].keys()))
    cell_list_len=len(cell_keys)
    datasets = [sim_data['compartment250k_mean1+var1'][f'cell_{i}'][()] for i in range(cell_list_len)]
    a = pd.DataFrame({f'cell_{i}': dataset for i, dataset in enumerate(datasets)})
    datasets = [sim_data['compartment250k_'+parm][f'cell_{i}'][()] for i in range(cell_list_len)]
    b = pd.DataFrame({f'cell_{i}': dataset for i, dataset in enumerate(datasets)})
    X.close()
    sim_data.close()
           
    # data1 = pd.read_csv(data1_path,index_col=0)
    # data2 = pd.read_csv(data2_path,index_col=0)

    # a=data1.iloc[:,:-5]
    # b=data2.iloc[:,:-5]
    #对于原始的高mean值的区间的bin的数据 使用改变后的mean值覆盖所有值，将该区间的Bin的方差置0
    #测试
    # threshold1 = 2
    # threshold2 = 3
    # red_points = (data2['Meanraw'].abs() > threshold1) & (data2['Meanraw'].abs() < threshold2)
    # b[red_points] = data2['Meannew'][red_points].values.reshape(-1, 1)
    #对于原始的高mean值的区间的bin的数据 使用改变后的mean值覆盖所有值，将该区间的Bin的方差置0
    #测试

    # c=data3.iloc[:,:-5]
    data=pd.merge(a,b,left_index=True,right_index=True)
    
    # data=pd.merge(data,c,left_index=True,right_index=True)
    label1=['A']*a.shape[1]
    label2=['B']*a.shape[1]
    # label3=['C']*a.shape[1]
    label=[]
    label.extend(label1)
    label.extend(label2)
    # label.extend(label3)
    dataset=load_sim_data(data,label,chr_index)
    # model= multi_input_cnn(input_shapes_list,label)
    model = Simple1DCNN(data,label)
    # model.summary()
    plot_model(model, to_file='model-chr1-chrX.png', show_shapes=True, show_layer_names=False, expand_nested=True)
    model = train(model,dataset,12,16,0.001)
    run_shap(model, data ,label,parm,chr_index)