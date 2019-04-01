import numpy as np
from tensorflow.keras.layers import Input, Dense, Lambda, Add, Layer
from tensorflow.keras import optimizers, models, regularizers
from tensorflow.keras import backend as K
import tensorflow as tf
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.models import load_model, Sequential

def resnet_for_dynamics(cf_trunc):

    num_epochs = 3000
    reg_param = 0.001
    
    features = np.transpose(cf_trunc)
    net_inputs = features[0:-1,:]
    net_outputs = features[1:,:]

    input_t = Input(shape=(np.shape(net_inputs)[1],))
     
    def input_wrapper(x):
        return x + input_t

    l1 = Dense(50, activation='tanh',kernel_regularizer=regularizers.l1(reg_param))(input_t)
    l1 = Dense(30, activation='tanh',kernel_regularizer=regularizers.l1(reg_param))(l1)
    l1 = Dense(10, activation='tanh',kernel_regularizer=regularizers.l1(reg_param))(l1)
    output = Dense(np.shape(net_inputs)[1], activation='linear')(l1)
    output_t = Lambda(lambda x: x * 1.0)(output)

    added = Lambda(lambda x: input_wrapper(x))(output_t)

    model = models.Model(inputs=input_t,outputs=added)
    my_adam = optimizers.Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.0, amsgrad=False)

    filepath = "best_weights_fcnn.h5"
    checkpoint = ModelCheckpoint(filepath, monitor='val_loss', verbose=1, save_best_only=True, mode='min',save_weights_only=True)
    callbacks_list = [checkpoint]

    model.compile(optimizer=my_adam,loss='mean_squared_error',metrics=[coeff_determination])
    history_callback = model.fit(net_inputs, net_outputs, epochs=num_epochs, batch_size=8, validation_split=0.33, callbacks=callbacks_list)

    model.load_weights(filepath)

    return model

def evaluate_rom_deployment(model,dataset,num_trunc,data_tsteps):
    # Make the initial condition from the first column of the dataset
    features = np.transpose(dataset)  
    input_state = np.copy(np.reshape(features[0,:],newshape=(1,2*num_trunc)))

    state_tracker = np.zeros(shape=(data_tsteps,np.shape(input_state)[1]),dtype='double')

    for t in range(data_tsteps):
        state_tracker[t,:] = input_state[:]
        output_state = model.predict(input_state)
        input_state = np.copy(output_state)

    return output_state, state_tracker

def coeff_determination(y_pred, y_true): #Order of function inputs is important here        
    SS_res =  K.sum(K.square( y_true-y_pred )) 
    SS_tot = K.sum(K.square( y_true - K.mean(y_true) ) )
    return ( 1 - SS_res/(SS_tot + K.epsilon()) )