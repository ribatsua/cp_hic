#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@File    :   main.py
@Time    :   2024/01/05 15:50:29
@Author  :   cp 
@Version :   1.0
@Desc    :   None
'''


import json
import shap
import h5py
import os
import numpy as np
from data_label_process import process_data, load_dataset, load_data_to_shap
from cnn_utils import  Simple1DCNN, train,multi_input_cnn
from tqdm import tqdm

np.random.seed(123)


def load_config(config_path):
    with open(config_path, 'r') as config_file:
        config = json.load(config_file)
    return config
def choice_data(X,y,back_name):
    data_list=[]
    for i in range(X.shape[0]):#数量应加1   Test_data.shape[0]+1
        # print("序号是：{} 标签是：{}".format(i-1,np.argmax(Test_lable[i-1])))
        if(np.argmax(y[i])==back_name):
            data_list.append(X[i])
    return np.stack(data_list)


def run_shap(model, X, y ,ori_label,data_dir,group):
    background=X[np.random.choice(X.shape[0], 100, replace=False)]
    #尝试分析两组数据比较，避开随机背景数据带来的结果变化
    # background=choice_data(X,ori_label,back_name=0)
    # ##################print('background shape:',background.shape)
    # a=X[0:343,:,:][np.random.choice(343, 50, replace=False)]
    # b=X[343:816,:,:][np.random.choice(473, 50, replace=False)]
    # c=X[816:1075,:,:][np.random.choice(259, 50, replace=False)]
    # background=np.concatenate((a,b,c),axis=0)
    # ex_data=choice_data(X,ori_label,back_name=1)
    # print('ex_data shape:',ex_data.shape)
    e = shap.DeepExplainer(model, background)
    print("e.expected_value:",e.expected_value)
    out_list = []
    num_samples = np.shape(X)[0]
    # num_samples = np.shape(ex_data)[0]
    # print("num_samples:",num_samples)
    for sample in tqdm(range(num_samples)):
        # shap
        shap_values = e.shap_values(X[sample : sample + 1])
        out_list.append(shap_values)
    shap_arr = np.squeeze(np.array(out_list))  #当模型为单输入时
    #当模型按染色体输入时，将各染色体shap值拼接在一起再处理 最终得到的shap_arr第一个维度是样本数量，第二个是输出类别数，第三个是特征数量
    # background = [array[np.random.choice(array.shape[0], 100, replace=False)] for array in X]
    # explainer = shap.GradientExplainer(model, background)
    # out_list = []
    # num_samples = X[0].shape[0]
    # # num_samples = np.shape(ex_data)[0]
    # # print("num_samples:",num_samples)
    # for sample in tqdm(range(num_samples)):
    #     # shap
    #     shap_values = explainer.shap_values([array[sample : sample + 1] for array in X])
    #     out_list.append(shap_values)
    # concatenated_third_level = [[np.concatenate(array_list, axis=1) for array_list in sublist] for sublist in out_list]
    # shap_arr = np.squeeze(np.array(concatenated_third_level))
    #当模型按染色体输入时，将各染色体shap值拼接在一起再处理 最终得到的shap_arr第一个维度是样本数量，第二个是输出类别数，第三个是特征数量

    save_path = os.path.join(data_dir,group)
    file_path= os.path.join(save_path,'compartment_250k_shap.h5')
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    with h5py.File(file_path, 'w') as hf:
        grp = hf.create_group(group)
        for i in range(shap_arr.shape[2]):
            # cell_data = shap_arr[i].T  
            bin_data=shap_arr[:,:,i] 
            # grp.create_dataset(f'cell_{i + 1}', data=cell_data)
            grp.create_dataset(f'bin_{i + 1}', data=bin_data)
        label_shap=grp.create_group('label')
        label_shap.create_dataset('cell_type', data=[l.encode('utf8') for l in y['cell_type']],
                            dtype=h5py.special_dtype(vlen=str))
        label_shap.create_dataset('cell_id', data=[l.encode('utf8') for l in y['cell_id']],
                            dtype=h5py.special_dtype(vlen=str))


    # for sample in tqdm(range(num_samples)):
    #     # shap
    #     shap_values = e.shap_values(ex_data[sample : sample + 1])
    #     out_list.append(shap_values)
    # shap_arr = np.squeeze(np.array(out_list))
    # save_path=data_dir+'/'+group+'_shap_values_150_remove_chrX_ncb.h5'
    
    # save_path = os.path.join(data_dir,group)
    # file_path= os.path.join(save_path,'compartemnt_250k_chr.h5')
    # if not os.path.exists(save_path):
    #     os.makedirs(save_path)
    # with h5py.File(file_path, 'w') as hf:
    #     grp = hf.create_group(group)
    #     for i in range(shap_arr.shape[2]):
    #         # cell_data = shap_arr[i].T  
    #         bin_data=shap_arr[:,:,i] 
    #         # grp.create_dataset(f'cell_{i + 1}', data=cell_data)
    #         grp.create_dataset(f'bin_{i + 1}', data=bin_data)
    #     label_shap=grp.create_group('label')
    #     label_shap.create_dataset('cell_type', data=[l.encode('utf8') for l in y['cell_type']],
    #                         dtype=h5py.special_dtype(vlen=str))
    #     label_shap.create_dataset('cell_id', data=[l.encode('utf8') for l in y['cell_id']],
    #                         dtype=h5py.special_dtype(vlen=str))
    
    hf.close()
    print("SHAP_value save to:",file_path)





def main():
    config_path='/home/script/python/config.json'
    config = load_config(config_path)
    data_dir = config['data_dir']
    mode = config['mode']
    group = config['group']
    epochs = config['epochs']
    batch_size = config['batch_size']
    lr = config['lr']
    chrom_list= config['chrom_list']
    smote = config['smote']
    print("data_dir:",data_dir,"\n","mode:",mode,"\n","group:",group,"\n","epochs:",epochs,"batch_size:",batch_size,"lr:",lr,"\n","chrom_list",chrom_list)

    data,label,chr_index =process_data(data_dir, mode, group,smote,chrom_list)
    # print(data)
    # dataset,input_shapes_list=load_dataset(data,label,chr_index)
    dataset=load_dataset(data,label,chr_index)
    model = Simple1DCNN(data,label)
    # model=train(model,dataset,epochs,batch_size,lr)
    # model = multi_input_cnn(input_shapes_list,label)

    model=train(model,dataset,epochs,batch_size,lr)
    ################特征冗余测试

    # model=load_model_test()
    ################特征冗余测试
    # X, y,ori_label,input_shapes_list=load_data_to_shap(data_dir, mode, group,model,chrom_list)
    X, y,ori_label=load_data_to_shap(data_dir, mode, group,model,chrom_list)

    run_shap(model,X, y ,ori_label,data_dir,group)


if __name__ == '__main__':
    main()
