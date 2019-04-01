import numpy as np
from tensorflow.keras.layers import Input, Dense, Lambda, Add, LSTM
from tensorflow.keras import optimizers, models, regularizers
from tensorflow.keras import backend as K
import tensorflow as tf
from tensorflow.keras.models import load_model, Sequential
from tensorflow.keras.callbacks import ModelCheckpoint

def lstm_for_dynamics(cf_trunc,seq_num=5):
    features = np.transpose(cf_trunc)
    states = np.copy(features[:,:]) #Rows are time, Columns are state values

    # Need to make batches of 10 input sequences and 1 output
    batch_size = np.shape(features)[0]-seq_num
    input_seq = np.zeros(shape=(batch_size,seq_num,np.shape(states)[1]))
    output_seq = np.zeros(shape=(batch_size,np.shape(states)[1]))

    for t in range(batch_size):
    	input_seq[t,:,:] = states[None,t:t+seq_num,:]
    	output_seq[t,:] = states[t+seq_num,:]
    
    # Model architecture
    model = Sequential()
    model.add(LSTM(32, return_sequences=True,
               input_shape=(seq_num, np.shape(states)[1])))  # returns a sequence of vectors of dimension 32
    model.add(LSTM(32))  # returns a sequence of vectors of dimension 32
    model.add(Dense(np.shape(states)[1], activation='linear'))

    num_epochs = 200

    # design network
    my_adam = optimizers.Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.0, amsgrad=False)

    filepath = "best_weights_lstm.h5"
    checkpoint = ModelCheckpoint(filepath, monitor='val_loss', verbose=1, save_best_only=True, mode='min',save_weights_only=True)
    callbacks_list = [checkpoint]
    
    # fit network
    model.compile(optimizer=my_adam,loss='mean_squared_error',metrics=[coeff_determination])
    train_history = model.fit(input_seq, output_seq, epochs=num_epochs, batch_size=8, validation_split=0.33, callbacks=callbacks_list)

    model.load_weights(filepath)

    return model


def evaluate_rom_deployment_lstm(model,dataset,num_trunc,data_tsteps,seq_num=5):

    # Make the initial condition from the first seq_num columns of the dataset
    features = np.transpose(dataset)  
    input_state = np.copy(np.reshape(features[0:seq_num,:],newshape=(seq_num,2*num_trunc)))

    state_tracker = np.zeros(shape=(1,data_tsteps,np.shape(features)[1]),dtype='double')
    state_tracker[0,0:seq_num,:] = input_state[0:seq_num,:]
    
    for t in range(seq_num,data_tsteps):
        lstm_input = state_tracker[:,t-seq_num:t,:]
        output_state = model.predict(lstm_input)
        state_tracker[0,t,:] = output_state[:]

    return output_state, state_tracker[0,:,:]

def coeff_determination(y_pred, y_true): #Order of function inputs is important here        
    SS_res =  K.sum(K.square( y_true-y_pred )) 
    SS_tot = K.sum(K.square( y_true - K.mean(y_true) ) )
    return ( 1 - SS_res/(SS_tot + K.epsilon()) )