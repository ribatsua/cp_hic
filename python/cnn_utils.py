#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@File    :   cnn_utils.py
@Time    :   2024/01/05 15:43:19
@Author  :   cp 
@Version :   1.0
@Desc    :   构建模型及训练部分代码
'''
# from keras.models import Sequential,load_model
# from keras.layers import Conv1D, MaxPooling1D, Flatten, Dense, Dropout
# from keras.optimizers import Adam
# from keras.callbacks import EarlyStopping
# Deep Neural Networks:
import tensorflow as tf; print('We\'re using TF-{}.'.format(tf.__version__))
# import keras; print('We\'re using Keras-{}.'.format(keras.__version__))
from keras.layers import (Input, Dense, Dropout, Flatten, BatchNormalization,
                                     Conv1D, Conv2D, MaxPooling1D, MaxPooling2D,
                                     LSTM, GRU, Embedding, Bidirectional, Concatenate,GlobalMaxPooling1D)
from keras.regularizers import (l1, l2, l1_l2)
from keras.optimizers import (RMSprop, Adam, SGD)
from keras.models import (Sequential, Model)

# Core:
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interp
import matplotlib.patches as patches

# Performance:
from sklearn.metrics import (confusion_matrix, classification_report, matthews_corrcoef, precision_score, roc_curve, auc)
from sklearn.model_selection import (StratifiedKFold, KFold, train_test_split)

#Utilities:
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.utils import to_categorical as labelEncoding
from keras.utils import plot_model  
 


def Simple1DCNN(data,label):
    bin_num=data.shape[0]
    # class_num=len(list(set(label['cell_type'])))
    class_num = 2
    model = Sequential()
    model.add(Conv1D(8,3,activation='relu',input_shape=(bin_num,1)))
    model.add(MaxPooling1D(2))
    model.add(Conv1D(16,3,activation='relu'))
    model.add(MaxPooling1D(2))
    model.add(Flatten())
    # model.add(Dense(16,activation='relu'))
    # model.add(Dropout(0.5))
    model.add(Dense(class_num, activation='softmax'))
    model.summary()
    return model

def train(model,dataset,epochs,batch_size,lr):
    Optimizer=Adam(learning_rate=lr)
    X_train_3d, y_train, X_test_3d ,y_test=dataset
    model.compile(loss='categorical_crossentropy', optimizer=Optimizer, 
                    metrics=['accuracy'])
    early_stopping = EarlyStopping(monitor='val_loss',  
                                mode='min',          
                                patience=5,          
                                verbose=1)

    model.fit(X_train_3d, y_train, epochs=epochs, batch_size=batch_size, validation_data=(X_test_3d, y_test), callbacks=[early_stopping])
    _, accuracy = model.evaluate(X_test_3d, y_test)
    print('Accuracy: %.2f%%' % (accuracy * 100))
    # model.save("test_model.h5")

    return model

# def load_model_test():
#     model=load_model("/home/test_model.h5")
#     #输出不作softmax尝试
#     # model.layers[-1].activation = None  # 设置输出层激活函数为 linear

#     # # 编译模型
#     # model.compile(optimizer=Adam(learning_rate=0.001), loss='categorical_crossentropy', metrics=['accuracy'])
#     return model

def Simple1DCNN_3d(data,label):
    bin_num=data.shape[1]
    class_num=len(list(set(label['cell_type'])))
    model = Sequential()
    model.add(Conv1D(8,3,activation='relu',input_shape=(bin_num,3)))
    model.add(MaxPooling1D(2))
    model.add(Conv1D(16,3,activation='relu'))
    model.add(MaxPooling1D(2))
    model.add(Flatten())
    # model.add(Dense(16,activation='relu'))
    # model.add(Dropout(0.5))
    model.add(Dense(class_num, activation='softmax'))
    model.summary()
    return model

def multi_input_cnn(input_shapes,label):
    inputs = []
    heads = [] 
    
    # class_num=2
                                                                                                           
    class_num=len(list(set(label['cell_type'])))

    # Create input layers for each head
    for input_shape in input_shapes:
        input_layer = Input(shape=input_shape)
        inputs.append(input_layer)
        x = Conv1D(filters=8, kernel_size=3, activation='relu')(input_layer)
        x = Conv1D(filters=8, kernel_size=3, activation='relu')(x)
        x = MaxPooling1D(2)(x)
        x = Conv1D(filters=16, kernel_size=3, activation='relu')(x)
        x = Conv1D(filters=16, kernel_size=3, activation='relu')(x)
        x = MaxPooling1D(2)(x)
        # x = GlobalMaxPooling1D()(x)
        head = Flatten()(x)
        heads.append(head)

    # Concatenate all heads
    merge = Concatenate()(heads)

    output = Dense(units=class_num, activation='softmax')(merge)
    return Model(inputs=inputs, outputs=output)

