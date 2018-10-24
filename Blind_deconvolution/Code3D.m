%Romit Maulik - CFD Lab
%Algorithms for training
%Levenberg Marquardt, Bayesian Regularization, Scaled Conjugate Gradient
%Extreme Learning Machine
%Trapezoidal filtering process.

toolbox_use = 2;%If 1- use nnstart else use ELM
toolbox_algo = 2;%If 1 use BR, 2 use LM, 3 use CG
hidden_neurons = 100;

noise_amp = 0.2;
noise_amp_test = 0.22;
sigma_val = 0.0;
sigma_val_test = 0.0;
stencil_type = 1;%1-For all neighbors, 2-For grid neighbors

%Importing coarse data files
raw_data = importdata('KHI_3D_Coarse.plt');
raw_data = raw_data.data;
test_data = importdata('KHI_3D_Coarse_Shifted.plt');
test_data = test_data.data;


%Storing data according fortran scheme i.e. omega(:,:)
%Defining vorticity
nx = 64;
u = zeros(nx+3,nx+3,nx+3);%2 is where the domain starts - 1 is ghost
v = zeros(nx+3,nx+3,nx+3);%2 is where the domain starts - 1 is ghost
w = zeros(nx+3,nx+3,nx+3);%2 is where the domain starts - 1 is ghost
%Importing into vort
iter = 1;
for k=2:nx+2
    for j=2:nx+2
        for i=2:nx+2
            u(i,j,k) = raw_data(iter,4);
            v(i,j,k) = raw_data(iter,5);
            w(i,j,k) = raw_data(iter,6);
            iter = iter + 1;
        end
    end
end

%Periodic boundary conditions
for k = 1:nx+3
for j = 1:nx+3
   u(1,j,k) = u(nx+2,j,k); 
   u(nx+3,j,k) = u(2,j,k);
   v(1,j,k) = v(nx+2,j,k); 
   v(nx+3,j,k) = v(2,j,k);
   w(1,j,k) = w(nx+2,j,k); 
   w(nx+3,j,k) = w(2,j,k);   
end
end

for k=1:nx+3
for i = 1:nx+3
   u(i,1,k) = u(i,nx+2,k); 
   u(i,nx+3,k) = u(i,2,k);
   v(i,1,k) = v(i,nx+2,k); 
   v(i,nx+3,k) = v(i,2,k);
   w(i,1,k) = w(i,nx+2,k); 
   w(i,nx+3,k) = w(i,2,k);   
end
end

for j=1:nx+3
for i=1:nx+3
   u(i,j,1) = u(i,j,nx+2); 
   u(i,j,nx+3) = u(i,j,2); 
   v(i,j,1) = v(i,j,nx+2); 
   v(i,j,nx+3) = v(i,j,2); 
   w(i,j,1) = w(i,j,nx+2); 
   w(i,j,nx+3) = w(i,j,2); 
end
end

%Normalizing u between -1 and 1
u = -1 + 2.*(u - min(u(:)))./(max(u(:)) - min(u(:)));
v = -1 + 2.*(v - min(v(:)))./(max(v(:)) - min(v(:)));
w = -1 + 2.*(w - min(w(:)))./(max(w(:)) - min(w(:)));

if sigma_val>0
uf = imgaussfilt3(u,sigma_val);
vf = imgaussfilt3(v,sigma_val);
wf = imgaussfilt3(w,sigma_val);
end
    
if sigma_val>0    
%Adding noise to uf,vf,wf
uf = uf + noise_amp.*randn(size(uf));
vf = vf + noise_amp.*randn(size(vf));
wf = wf + noise_amp.*randn(size(wf));
else
uf = u + noise_amp.*randn(size(u));
vf = v + noise_amp.*randn(size(v));
wf = w + noise_amp.*randn(size(w));
end


%Training for u
if stencil_type==1
    inputs = zeros(27,64^3);
else
    inputs = zeros(7,64^3);
end

outputs =zeros(1,64^3);


%Prepare datasets for NN training
iter = 1;

for k=2:nx+2
for j=2:nx+2
    for i = 2:nx+2
        
    if stencil_type == 1
        inputs(1,iter) = uf(i,j,k);
        inputs(2,iter) = uf(i,j,k+1);
        inputs(3,iter) = uf(i,j,k-1);
        inputs(4,iter) = uf(i,j+1,k);
        inputs(5,iter) = uf(i,j+1,k+1);
        inputs(6,iter) = uf(i,j+1,k-1);
        inputs(7,iter) = uf(i,j-1,k);
        inputs(8,iter) = uf(i,j-1,k+1);
        inputs(9,iter) = uf(i,j-1,k-1);        
        
        inputs(10,iter) = uf(i+1,j,k);
        inputs(11,iter) = uf(i+1,j,k+1);
        inputs(12,iter) = uf(i+1,j,k-1);
        inputs(13,iter) = uf(i+1,j+1,k);
        inputs(14,iter) = uf(i+1,j+1,k+1);
        inputs(15,iter) = uf(i+1,j+1,k-1);
        inputs(16,iter) = uf(i+1,j-1,k);
        inputs(17,iter) = uf(i+1,j-1,k+1);
        inputs(18,iter) = uf(i+1,j-1,k-1);        
        
        inputs(19,iter) = uf(i-1,j,k);
        inputs(20,iter) = uf(i-1,j,k+1);
        inputs(21,iter) = uf(i-1,j,k-1);
        inputs(22,iter) = uf(i-1,j+1,k);
        inputs(23,iter) = uf(i-1,j+1,k+1);
        inputs(24,iter) = uf(i-1,j+1,k-1);
        inputs(25,iter) = uf(i-1,j-1,k);
        inputs(26,iter) = uf(i-1,j-1,k+1);
        inputs(27,iter) = uf(i-1,j-1,k-1);        
    else      
        inputs(1,iter) = uf(i,j,k);
        inputs(2,iter) = uf(i,j,k+1);
        inputs(3,iter) = uf(i,j,k-1);
        inputs(4,iter) = uf(i,j+1,k);
        inputs(5,iter) = uf(i,j-1,k);
        inputs(6,iter) = uf(i+1,j,k);
        inputs(7,iter) = uf(i-1,j,k);                      
    end
               
        outputs(1,iter) = u(i,j,k);
        iter = iter + 1;
    end
end
end

x = inputs;
t = outputs;


%Preparing test data for ELM use
%Storing data according fortran scheme i.e. omega(:,:)
%Defining vorticity
nx = 64;
u_test = zeros(nx+3,nx+3,nx+3);%2 is where the domain starts - 1 is ghost
v_test = zeros(nx+3,nx+3,nx+3);%2 is where the domain starts - 1 is ghost
w_test = zeros(nx+3,nx+3,nx+3);%2 is where the domain starts - 1 is ghost

%Importing into vort
iter = 1;
for k=2:nx+2
    for j=2:nx+2
        for i=2:nx+2
            u_test(i,j,k) = test_data(iter,4);
            v_test(i,j,k) = test_data(iter,5);
            w_test(i,j,k) = test_data(iter,6);
            iter = iter + 1;
        end
    end
end

%Periodic boundary conditions
for k = 1:nx+3
for j = 1:nx+3
   u_test(1,j,k) = u_test(nx+2,j,k); 
   u_test(nx+3,j,k) = u_test(2,j,k);
   v_test(1,j,k) = v_test(nx+2,j,k); 
   v_test(nx+3,j,k) = v_test(2,j,k);
   w_test(1,j,k) = w_test(nx+2,j,k); 
   w_test(nx+3,j,k) = w_test(2,j,k);   
end
end

for k=1:nx+3
for i = 1:nx+3
   u_test(i,1,k) = u_test(i,nx+2,k); 
   u_test(i,nx+3,k) = u_test(i,2,k);
   v_test(i,1,k) = v_test(i,nx+2,k); 
   v_test(i,nx+3,k) = v_test(i,2,k);
   w_test(i,1,k) = w_test(i,nx+2,k); 
   w_test(i,nx+3,k) = w_test(i,2,k);   
end
end

for j=1:nx+3
for i=1:nx+3
   u_test(i,j,1) = u_test(i,j,nx+2); 
   u_test(i,j,nx+3) = u_test(i,j,2); 
   v_test(i,j,1) = v_test(i,j,nx+2); 
   v_test(i,j,nx+3) = v_test(i,j,2); 
   w_test(i,j,1) = w_test(i,j,nx+2); 
   w_test(i,j,nx+3) = w_test(i,j,2); 
end
end

%Normalizing u between -1 and 1
u_test = -1 + 2.*(u_test - min(u_test(:)))./(max(u_test(:)) - min(u_test(:)));
v_test = -1 + 2.*(v_test - min(v_test(:)))./(max(v_test(:)) - min(v_test(:)));
w_test = -1 + 2.*(w_test - min(w_test(:)))./(max(w_test(:)) - min(w_test(:)));

if sigma_val>0
uf_test = imgaussfilt3(u_test,sigma_val_test);
vf_test = imgaussfilt3(v_test,sigma_val_test);
wf_test = imgaussfilt3(w_test,sigma_val_test);
end
    
if sigma_val>0    
%Adding noise to uf,vf,wf
uf_test = uf_test + noise_amp_test.*randn(size(uf_test));
vf_test = vf_test + noise_amp_test.*randn(size(vf_test));
wf_test = wf_test + noise_amp_test.*randn(size(wf_test));
else
uf_test = u_test + noise_amp_test.*randn(size(u_test));
vf_test = v_test + noise_amp_test.*randn(size(v_test));
wf_test = w_test + noise_amp_test.*randn(size(w_test));
end



%Prepare datasets for testing
iter = 1;

for k=2:nx+2
for j=2:nx+2
    for i = 2:nx+2
        
    if stencil_type == 1
        inputs_test(1,iter) = uf_test(i,j,k);
        inputs_test(2,iter) = uf_test(i,j,k+1);
        inputs_test(3,iter) = uf_test(i,j,k-1);
        inputs_test(4,iter) = uf_test(i,j+1,k);
        inputs_test(5,iter) = uf_test(i,j+1,k+1);
        inputs_test(6,iter) = uf_test(i,j+1,k-1);
        inputs_test(7,iter) = uf_test(i,j-1,k);
        inputs_test(8,iter) = uf_test(i,j-1,k+1);
        inputs_test(9,iter) = uf_test(i,j-1,k-1);        
        
        inputs_test(10,iter) = uf_test(i+1,j,k);
        inputs_test(11,iter) = uf_test(i+1,j,k+1);
        inputs_test(12,iter) = uf_test(i+1,j,k-1);
        inputs_test(13,iter) = uf_test(i+1,j+1,k);
        inputs_test(14,iter) = uf_test(i+1,j+1,k+1);
        inputs_test(15,iter) = uf_test(i+1,j+1,k-1);
        inputs_test(16,iter) = uf_test(i+1,j-1,k);
        inputs_test(17,iter) = uf_test(i+1,j-1,k+1);
        inputs_test(18,iter) = uf_test(i+1,j-1,k-1);        
        
        inputs_test(19,iter) = uf_test(i-1,j,k);
        inputs_test(20,iter) = uf_test(i-1,j,k+1);
        inputs_test(21,iter) = uf_test(i-1,j,k-1);
        inputs_test(22,iter) = uf_test(i-1,j+1,k);
        inputs_test(23,iter) = uf_test(i-1,j+1,k+1);
        inputs_test(24,iter) = uf_test(i-1,j+1,k-1);
        inputs_test(25,iter) = uf_test(i-1,j-1,k);
        inputs_test(26,iter) = uf_test(i-1,j-1,k+1);
        inputs_test(27,iter) = uf_test(i-1,j-1,k-1);        
    else      
        inputs_test(1,iter) = uf_test(i,j,k);
        inputs_test(2,iter) = uf_test(i,j,k+1);
        inputs_test(3,iter) = uf_test(i,j,k-1);
        inputs_test(4,iter) = uf_test(i,j+1,k);
        inputs_test(5,iter) = uf_test(i,j-1,k);
        inputs_test(6,iter) = uf_test(i+1,j,k);
        inputs_test(7,iter) = uf_test(i-1,j,k);                      
    end
               
        iter = iter + 1;
    end
end
end

if toolbox_use==1

% Choose a Training Function
% For a list of all training functions type: help nntrain
% 'trainlm' is usually fastest.
% 'trainbr' takes longer but may be better for challenging problems.
% 'trainscg' uses less memory. Suitable in low memory situations.

if toolbox_algo == 1
    trainFcn = 'trainbr';  % Bayesian Regularization backpropagation.
elseif toolbox_algo == 2
    trainFcn = 'trainlm';  %  Levenberg Marquardt backpropagation.
elseif toolbox_algo == 3
    trainFcn = 'trainscg';  % Conjugate Gradient backpropagation.
end
 
% Create a Fitting Network
hiddenLayerSize = hidden_neurons;
net = fitnet(hiddenLayerSize,trainFcn);

% Choose Input and Output Pre/Post-Processing Functions
% For a list of all processing functions type: help nnprocess
net.input.processFcns = {'removeconstantrows','mapminmax'};
net.output.processFcns = {'removeconstantrows','mapminmax'};

% Setup Division of Data for Training, Validation, Testing
% For a list of all data division functions type: help nndivide
net.divideFcn = 'dividerand';  % Divide data randomly
net.divideMode = 'sample';  % Divide up every sample
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;

% Choose a Performance Function
% For a list of all performance functions type: help nnperformance
net.performFcn = 'mse';  % Mean Squared Error

% Choose Plot Functions
% For a list of all plot functions type: help nnplot
net.plotFcns = {'plotperform','plottrainstate','ploterrhist', ...
    'plotregression', 'plotfit'};

% Train the Network
[net,tr] = train(net,x,t);

% Test the Network
y = net(x);
e = gsubtract(t,y);
performance = perform(net,t,y);

% Recalculate Training, Validation and Test Performance
trainTargets = t .* tr.trainMask{1};
%valTargets = t .* tr.valMask{1};%not needed for BR
testTargets = t .* tr.testMask{1};
trainPerformance = perform(net,trainTargets,y);
%valPerformance = perform(net,valTargets,y);%not needed for BR
testPerformance = perform(net,testTargets,y);

% View the Network
%view(net)

% Plots
% Uncomment these lines to enable various plots.
%figure, plotperform(tr)
%figure, plottrainstate(tr)
%figure, ploterrhist(e)
%figure, plotregression(t,y)
%figure, plotfit(net,x,t)

%Generate the function
if toolbox_algo==1
genFunction(net,'BR_Trained_Func_u');
NN_prediction_u = BR_Trained_Func_u(x);
elseif toolbox_algo==2
genFunction(net,'LM_Trained_Func_u');
NN_prediction_u = LM_Trained_Func_u(x);
elseif toolbox_algo==3
genFunction(net,'CG_Trained_Func_u');
NN_prediction_u = CG_Trained_Func_u(x);
end


else %Extreme Learning Machine Huang et al. Neurocomputing 70 - 2006
    
W1 = -1 + (1+1).*rand(hidden_neurons,size(inputs,1));
B1 = -1 + (1+1).*rand(hidden_neurons,1);

H = (W1*inputs);

Htest=(W1*inputs_test);

for j = 1:size(H,2)
    for i = 1:size(H,1)
        H(i,j) = H(i,j)+B1(i,1);
        H(i,j)=tansig(H(i,j));
    end
end

for j = 1:size(Htest,2)
    for i = 1:size(Htest,1)
        Htest(i,j) = Htest(i,j)+B1(i,1);
        Htest(i,j)=tansig(Htest(i,j));
    end
end

H = H';
Htest = Htest';
   
%calculating beta matrix
B = pinv(H)*t';

%Checking training
NN_prediction_u = Htest*B;

NN_prediction_u = NN_prediction_u';

end

%Training for v
if stencil_type==1
    inputs = zeros(27,64^3);
else
    inputs = zeros(7,64^3);
end

outputs =zeros(1,64^3);

%Prepare datasets for NN training
iter = 1;

for k=2:nx+2
for j=2:nx+2
    for i = 2:nx+2
        
        if stencil_type == 1
        inputs(1,iter) = vf(i,j,k);
        inputs(2,iter) = vf(i,j,k+1);
        inputs(3,iter) = vf(i,j,k-1);
        inputs(4,iter) = vf(i,j+1,k);
        inputs(5,iter) = vf(i,j+1,k+1);
        inputs(6,iter) = vf(i,j+1,k-1);
        inputs(7,iter) = vf(i,j-1,k);
        inputs(8,iter) = vf(i,j-1,k+1);
        inputs(9,iter) = vf(i,j-1,k-1);        
        
        inputs(10,iter) = vf(i+1,j,k);
        inputs(11,iter) = vf(i+1,j,k+1);
        inputs(12,iter) = vf(i+1,j,k-1);
        inputs(13,iter) = vf(i+1,j+1,k);
        inputs(14,iter) = vf(i+1,j+1,k+1);
        inputs(15,iter) = vf(i+1,j+1,k-1);
        inputs(16,iter) = vf(i+1,j-1,k);
        inputs(17,iter) = vf(i+1,j-1,k+1);
        inputs(18,iter) = vf(i+1,j-1,k-1);        
        
        inputs(19,iter) = vf(i-1,j,k);
        inputs(20,iter) = vf(i-1,j,k+1);
        inputs(21,iter) = vf(i-1,j,k-1);
        inputs(22,iter) = vf(i-1,j+1,k);
        inputs(23,iter) = vf(i-1,j+1,k+1);
        inputs(24,iter) = vf(i-1,j+1,k-1);
        inputs(25,iter) = vf(i-1,j-1,k);
        inputs(26,iter) = vf(i-1,j-1,k+1);
        inputs(27,iter) = vf(i-1,j-1,k-1);
        
        else
            
        inputs(1,iter) = vf(i,j,k);
        inputs(2,iter) = vf(i,j,k+1);
        inputs(3,iter) = vf(i,j,k-1);
        inputs(4,iter) = vf(i,j+1,k);
        inputs(5,iter) = vf(i,j-1,k);
        inputs(6,iter) = vf(i+1,j,k);
        inputs(7,iter) = vf(i-1,j,k);            
            
        end
               
        outputs(1,iter) = v(i,j,k);
        iter = iter + 1;
    end
end
end

x = inputs;
t = outputs;


%Prepare datasets for testing
iter = 1;

for k=2:nx+2
for j=2:nx+2
    for i = 2:nx+2
        
    if stencil_type == 1
        inputs_test(1,iter) = vf_test(i,j,k);
        inputs_test(2,iter) = vf_test(i,j,k+1);
        inputs_test(3,iter) = vf_test(i,j,k-1);
        inputs_test(4,iter) = vf_test(i,j+1,k);
        inputs_test(5,iter) = vf_test(i,j+1,k+1);
        inputs_test(6,iter) = vf_test(i,j+1,k-1);
        inputs_test(7,iter) = vf_test(i,j-1,k);
        inputs_test(8,iter) = vf_test(i,j-1,k+1);
        inputs_test(9,iter) = vf_test(i,j-1,k-1);        
        
        inputs_test(10,iter) = vf_test(i+1,j,k);
        inputs_test(11,iter) = vf_test(i+1,j,k+1);
        inputs_test(12,iter) = vf_test(i+1,j,k-1);
        inputs_test(13,iter) = vf_test(i+1,j+1,k);
        inputs_test(14,iter) = vf_test(i+1,j+1,k+1);
        inputs_test(15,iter) = vf_test(i+1,j+1,k-1);
        inputs_test(16,iter) = vf_test(i+1,j-1,k);
        inputs_test(17,iter) = vf_test(i+1,j-1,k+1);
        inputs_test(18,iter) = vf_test(i+1,j-1,k-1);        
        
        inputs_test(19,iter) = vf_test(i-1,j,k);
        inputs_test(20,iter) = vf_test(i-1,j,k+1);
        inputs_test(21,iter) = vf_test(i-1,j,k-1);
        inputs_test(22,iter) = vf_test(i-1,j+1,k);
        inputs_test(23,iter) = vf_test(i-1,j+1,k+1);
        inputs_test(24,iter) = vf_test(i-1,j+1,k-1);
        inputs_test(25,iter) = vf_test(i-1,j-1,k);
        inputs_test(26,iter) = vf_test(i-1,j-1,k+1);
        inputs_test(27,iter) = vf_test(i-1,j-1,k-1);        
    else      
        inputs_test(1,iter) = vf_test(i,j,k);
        inputs_test(2,iter) = vf_test(i,j,k+1);
        inputs_test(3,iter) = vf_test(i,j,k-1);
        inputs_test(4,iter) = vf_test(i,j+1,k);
        inputs_test(5,iter) = vf_test(i,j-1,k);
        inputs_test(6,iter) = vf_test(i+1,j,k);
        inputs_test(7,iter) = vf_test(i-1,j,k);                      
    end
               
        iter = iter + 1;
    end
end
end

if toolbox_use==1

% Choose a Training Function
% For a list of all training functions type: help nntrain
% 'trainlm' is usually fastest.
% 'trainbr' takes longer but may be better for challenging problems.
% 'trainscg' uses less memory. Suitable in low memory situations.

if toolbox_algo == 1
    trainFcn = 'trainbr';  % Bayesian Regularization backpropagation.
elseif toolbox_algo == 2
    trainFcn = 'trainlm';  %  Levenberg Marquardt backpropagation.
elseif toolbox_algo == 3
    trainFcn = 'trainscg';  % Conjugate Gradient backpropagation.
end
 
% Create a Fitting Network
hiddenLayerSize = hidden_neurons;
net = fitnet(hiddenLayerSize,trainFcn);

% Choose Input and Output Pre/Post-Processing Functions
% For a list of all processing functions type: help nnprocess
net.input.processFcns = {'removeconstantrows','mapminmax'};
net.output.processFcns = {'removeconstantrows','mapminmax'};

% Setup Division of Data for Training, Validation, Testing
% For a list of all data division functions type: help nndivide
net.divideFcn = 'dividerand';  % Divide data randomly
net.divideMode = 'sample';  % Divide up every sample
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;

% Choose a Performance Function
% For a list of all performance functions type: help nnperformance
net.performFcn = 'mse';  % Mean Squared Error

% Choose Plot Functions
% For a list of all plot functions type: help nnplot
net.plotFcns = {'plotperform','plottrainstate','ploterrhist', ...
    'plotregression', 'plotfit'};

% Train the Network
[net,tr] = train(net,x,t);

% Test the Network
y = net(x);
e = gsubtract(t,y);
performance = perform(net,t,y);

% Recalculate Training, Validation and Test Performance
trainTargets = t .* tr.trainMask{1};
%valTargets = t .* tr.valMask{1};%not needed for BR
testTargets = t .* tr.testMask{1};
trainPerformance = perform(net,trainTargets,y);
%valPerformance = perform(net,valTargets,y);%not needed for BR
testPerformance = perform(net,testTargets,y);

% View the Network
%view(net)

% Plots
% Uncomment these lines to enable various plots.
%figure, plotperform(tr)
%figure, plottrainstate(tr)
%figure, ploterrhist(e)
%figure, plotregression(t,y)
%figure, plotfit(net,x,t)

%Generate the function
if toolbox_algo==1
genFunction(net,'BR_Trained_Func_v');
NN_prediction_v = BR_Trained_Func_v(inputs_test);
elseif toolbox_algo==2
genFunction(net,'LM_Trained_Func_v');
NN_prediction_v = LM_Trained_Func_v(inputs_test);
elseif toolbox_algo==3
genFunction(net,'CG_Trained_Func_v');
NN_prediction_v = CG_Trained_Func_v(inputs_test);
end


else %Extreme Learning Machine Huang et al. Neurocomputing 70 - 2006
    
W1 = -1 + (1+1).*rand(hidden_neurons,size(inputs,1));
B1 = -1 + (1+1).*rand(hidden_neurons,1);

H = (W1*inputs);
Htest = (W1*inputs_test);

for j = 1:size(H,2)
    for i = 1:size(H,1)
        H(i,j) = H(i,j)+B1(i,1);
        H(i,j)=tansig(H(i,j));
    end
end

for j = 1:size(Htest,2)
    for i = 1:size(Htest,1)
        Htest(i,j) = Htest(i,j)+B1(i,1);
        Htest(i,j)=tansig(Htest(i,j));
    end
end


H = H';
Htest = Htest';
   
%calculating beta matrix
B = pinv(H)*t';

%Checking training
NN_prediction_v = Htest*B;

NN_prediction_v = NN_prediction_v';

end



%Training for w
if stencil_type==1
    inputs = zeros(27,64^3);
else
    inputs = zeros(7,64^3);
end

%Prepare datasets for NN training
iter = 1;

for k=2:nx+2
for j=2:nx+2
    for i = 2:nx+2
        
        if stencil_type == 1
        inputs(1,iter) = wf(i,j,k);
        inputs(2,iter) = wf(i,j,k+1);
        inputs(3,iter) = wf(i,j,k-1);
        inputs(4,iter) = wf(i,j+1,k);
        inputs(5,iter) = wf(i,j+1,k+1);
        inputs(6,iter) = wf(i,j+1,k-1);
        inputs(7,iter) = wf(i,j-1,k);
        inputs(8,iter) = wf(i,j-1,k+1);
        inputs(9,iter) = wf(i,j-1,k-1);        
        
        inputs(10,iter) = wf(i+1,j,k);
        inputs(11,iter) = wf(i+1,j,k+1);
        inputs(12,iter) = wf(i+1,j,k-1);
        inputs(13,iter) = wf(i+1,j+1,k);
        inputs(14,iter) = wf(i+1,j+1,k+1);
        inputs(15,iter) = wf(i+1,j+1,k-1);
        inputs(16,iter) = wf(i+1,j-1,k);
        inputs(17,iter) = wf(i+1,j-1,k+1);
        inputs(18,iter) = wf(i+1,j-1,k-1);        
        
        inputs(19,iter) = wf(i-1,j,k);
        inputs(20,iter) = wf(i-1,j,k+1);
        inputs(21,iter) = wf(i-1,j,k-1);
        inputs(22,iter) = wf(i-1,j+1,k);
        inputs(23,iter) = wf(i-1,j+1,k+1);
        inputs(24,iter) = wf(i-1,j+1,k-1);
        inputs(25,iter) = wf(i-1,j-1,k);
        inputs(26,iter) = wf(i-1,j-1,k+1);
        inputs(27,iter) = wf(i-1,j-1,k-1);
        
        
        else
            
        inputs(1,iter) = wf(i,j,k);
        inputs(2,iter) = wf(i,j,k+1);
        inputs(3,iter) = wf(i,j,k-1);
        inputs(4,iter) = wf(i,j+1,k);
        inputs(5,iter) = wf(i,j-1,k);
        inputs(6,iter) = wf(i+1,j,k);
        inputs(7,iter) = wf(i-1,j,k);                      
        
        end
               
        outputs(1,iter) = w(i,j,k);
        iter = iter + 1;
    end
end
end

x = inputs;
t = outputs;

%Prepare datasets for testing
iter = 1;

for k=2:nx+2
for j=2:nx+2
    for i = 2:nx+2
        
    if stencil_type == 1
        inputs_test(1,iter) = wf_test(i,j,k);
        inputs_test(2,iter) = wf_test(i,j,k+1);
        inputs_test(3,iter) = wf_test(i,j,k-1);
        inputs_test(4,iter) = wf_test(i,j+1,k);
        inputs_test(5,iter) = wf_test(i,j+1,k+1);
        inputs_test(6,iter) = wf_test(i,j+1,k-1);
        inputs_test(7,iter) = wf_test(i,j-1,k);
        inputs_test(8,iter) = wf_test(i,j-1,k+1);
        inputs_test(9,iter) = wf_test(i,j-1,k-1);        
        
        inputs_test(10,iter) = wf_test(i+1,j,k);
        inputs_test(11,iter) = wf_test(i+1,j,k+1);
        inputs_test(12,iter) = wf_test(i+1,j,k-1);
        inputs_test(13,iter) = wf_test(i+1,j+1,k);
        inputs_test(14,iter) = wf_test(i+1,j+1,k+1);
        inputs_test(15,iter) = wf_test(i+1,j+1,k-1);
        inputs_test(16,iter) = wf_test(i+1,j-1,k);
        inputs_test(17,iter) = wf_test(i+1,j-1,k+1);
        inputs_test(18,iter) = wf_test(i+1,j-1,k-1);        
        
        inputs_test(19,iter) = wf_test(i-1,j,k);
        inputs_test(20,iter) = wf_test(i-1,j,k+1);
        inputs_test(21,iter) = wf_test(i-1,j,k-1);
        inputs_test(22,iter) = wf_test(i-1,j+1,k);
        inputs_test(23,iter) = wf_test(i-1,j+1,k+1);
        inputs_test(24,iter) = wf_test(i-1,j+1,k-1);
        inputs_test(25,iter) = wf_test(i-1,j-1,k);
        inputs_test(26,iter) = wf_test(i-1,j-1,k+1);
        inputs_test(27,iter) = wf_test(i-1,j-1,k-1);        
    else      
        inputs_test(1,iter) = wf_test(i,j,k);
        inputs_test(2,iter) = wf_test(i,j,k+1);
        inputs_test(3,iter) = wf_test(i,j,k-1);
        inputs_test(4,iter) = wf_test(i,j+1,k);
        inputs_test(5,iter) = wf_test(i,j-1,k);
        inputs_test(6,iter) = wf_test(i+1,j,k);
        inputs_test(7,iter) = wf_test(i-1,j,k);                      
    end
        iter = iter + 1;
    end
end
end

if toolbox_use==1

% Choose a Training Function
% For a list of all training functions type: help nntrain
% 'trainlm' is usually fastest.
% 'trainbr' takes longer but may be better for challenging problems.
% 'trainscg' uses less memory. Suitable in low memory situations.

if toolbox_algo == 1
    trainFcn = 'trainbr';  % Bayesian Regularization backpropagation.
elseif toolbox_algo == 2
    trainFcn = 'trainlm';  %  Levenberg Marquardt backpropagation.
elseif toolbox_algo == 3
    trainFcn = 'trainscg';  % Conjugate Gradient backpropagation.
end
 
% Create a Fitting Network
hiddenLayerSize = hidden_neurons;
net = fitnet(hiddenLayerSize,trainFcn);

% Choose Input and Output Pre/Post-Processing Functions
% For a list of all processing functions type: help nnprocess
net.input.processFcns = {'removeconstantrows','mapminmax'};
net.output.processFcns = {'removeconstantrows','mapminmax'};

% Setup Division of Data for Training, Validation, Testing
% For a list of all data division functions type: help nndivide
net.divideFcn = 'dividerand';  % Divide data randomly
net.divideMode = 'sample';  % Divide up every sample
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;

% Choose a Performance Function
% For a list of all performance functions type: help nnperformance
net.performFcn = 'mse';  % Mean Squared Error

% Choose Plot Functions
% For a list of all plot functions type: help nnplot
net.plotFcns = {'plotperform','plottrainstate','ploterrhist', ...
    'plotregression', 'plotfit'};

% Train the Network
[net,tr] = train(net,x,t);

% Test the Network
y = net(x);
e = gsubtract(t,y);
performance = perform(net,t,y);

% Recalculate Training, Validation and Test Performance
trainTargets = t .* tr.trainMask{1};
%valTargets = t .* tr.valMask{1};%not needed for BR
testTargets = t .* tr.testMask{1};
trainPerformance = perform(net,trainTargets,y);
%valPerformance = perform(net,valTargets,y);%not needed for BR
testPerformance = perform(net,testTargets,y);

% View the Network
%view(net)

% Plots
% Uncomment these lines to enable various plots.
%figure, plotperform(tr)
%figure, plottrainstate(tr)
%figure, ploterrhist(e)
%figure, plotregression(t,y)
%figure, plotfit(net,x,t)

%Generate the function
if toolbox_algo==1
genFunction(net,'BR_Trained_Func_w');
NN_prediction_w = BR_Trained_Func_w(inputs_test);
elseif toolbox_algo==2
genFunction(net,'LM_Trained_Func_w');
NN_prediction_w = LM_Trained_Func_w(inputs_test);
elseif toolbox_algo==3
genFunction(net,'CG_Trained_Func_w');
NN_prediction_w = CG_Trained_Func_w(inputs_test);
end


else %Extreme Learning Machine Huang et al. Neurocomputing 70 - 2006
    
W1 = -1 + (1+1).*rand(hidden_neurons,size(inputs,1));
B1 = -1 + (1+1).*rand(hidden_neurons,1);

H = (W1*inputs);
Htest = (W1*inputs_test);

for j = 1:size(H,2)
    for i = 1:size(H,1)
        H(i,j) = H(i,j)+B1(i,1);
        H(i,j)=tansig(H(i,j));
    end
end

for j = 1:size(Htest,2)
    for i = 1:size(Htest,1)
        Htest(i,j) = Htest(i,j)+B1(i,1);
        Htest(i,j)=tansig(Htest(i,j));
    end
end

H = H';
Htest = Htest';
   
%calculating beta matrix
B = pinv(H)*t';

%Checking training
NN_prediction_w = Htest*B;

NN_prediction_w = NN_prediction_w';

end

%Exporting into nn arrays
unn = zeros(nx+3,nx+3,nx+3);%2 is where the domain starts - 1 is ghost
vnn = zeros(nx+3,nx+3,nx+3);%2 is where the domain starts - 1 is ghost
wnn = zeros(nx+3,nx+3,nx+3);%2 is where the domain starts - 1 is ghost
iter = 1;
for k=2:nx+2%Must be in k-i-j order
for j=2:nx+2
    for i = 2:nx+2      
        unn(i,j,k) = NN_prediction_u(1,iter);
        vnn(i,j,k) = NN_prediction_v(1,iter);
        wnn(i,j,k) = NN_prediction_w(1,iter);
        iter = iter + 1;
    end
end
end

%Outputting NN_predicted data
iter = 1;
for k=2:nx+2
for j=2:nx+2
    for i=2:nx+2
        raw_data(iter,4) = unn(i,j,k);
        raw_data(iter,5) = vnn(i,j,k);
        raw_data(iter,6) = wnn(i,j,k);
        iter = iter + 1;
    end
end
end

raw_data_op_new = raw_data';

fileID=fopen('nn_velocity.plt','w');
fprintf(fileID,' variables = "x", "y", "z", "u","v","w"\n');
fprintf(fileID,'zone T=Zone_4 i=          65  j=          65  k=          65  f=point\n');

formatSpec = '%16.14f         %16.14f         %16.14f         %16.14f         %16.14f         %16.14f\n';
fprintf(fileID,formatSpec,raw_data_op_new);

fclose(fileID);

%Outputting noisy filtered data
iter = 1;
for k=2:nx+2
for j=2:nx+2
    for i=2:nx+2
        raw_data(iter,4) = uf_test(i,j,k);
        raw_data(iter,5) = vf_test(i,j,k);
        raw_data(iter,6) = wf_test(i,j,k);
        iter = iter + 1;
    end
end
end

raw_data_op_new = raw_data';

fileID=fopen('velocity_f_noise.plt','w');
fprintf(fileID,' variables = "x", "y", "z", "u","v","w"\n');
fprintf(fileID,'zone T=Zone_4 i=          65  j=          65  k=          65  f=point\n');

formatSpec = '%16.14f         %16.14f         %16.14f         %16.14f         %16.14f         %16.14f\n';
fprintf(fileID,formatSpec,raw_data_op_new);

fclose(fileID);


%Outputting normalized input data
iter = 1;
for k=2:nx+2
for j=2:nx+2
    for i=2:nx+2
        raw_data(iter,4) = u_test(i,j,k);
        raw_data(iter,5) = v_test(i,j,k);
        raw_data(iter,6) = w_test(i,j,k);
        iter = iter + 1;
    end
end
end

raw_data_op_new = raw_data';

fileID=fopen('velocity_coarse.plt','w');
fprintf(fileID,' variables = "x", "y", "z", "u","v","w"\n');
fprintf(fileID,'zone T=Zone_4 i=          65  j=          65  k=          65  f=point\n');

formatSpec = '%16.14f         %16.14f         %16.14f         %16.14f         %16.14f         %16.14f\n';
fprintf(fileID,formatSpec,raw_data_op_new);

fclose(fileID);
