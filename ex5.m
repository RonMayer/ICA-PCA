%% Q1

% load files
%mixed files
clear;
clc;

Mfiles = {'mix1.wav' 'mix2.wav' 'mix3.wav'};
for i=1:length(Mfiles)
    [My(:,i),Mfs(:,i)] = audioread(Mfiles{i});
end
figure()
for i=1:size(My,2)
    subplot(3,1,i);
    plot(My(:,i)');
    title(['audio # ' num2str(i)]);
    hold on
end
sgtitle('mixed files');

% source files

Sfiles = {'source1.wav' 'source2.wav' 'source3.wav'};
for i=1:length(Sfiles)
    [Sy(:,i),Sfs(:,i)] = audioread(Sfiles{i});
end
figure()
for i=1:size(Sy,2)
    subplot(3,1,i);
    plot(Sy(:,i)');
    title(['audio # ' num2str(i)]);
    hold on
end
sgtitle('source files');    

% infomaxICA

r=size(Sy,2); % #of independent components
eta=0.3; %learning rate
iterMax=10; %max iterations
conLim=1e-5; %convergence criteria

%initialize W
normRows = @(X) bsxfun(@rdivide,X,sqrt(sum(X.^2,2)));
W = normRows(rand(r,size(My,2))); % Random initial weights
delta=1;
k=0;

while delta > conLim && k < iterMax
    k=k+1;
    Wlast=W;
    sample = randperm(length(My));
    for t=1:length(My)
        x = My(sample(t),:)';
        y=W*x; 
        dW=(eta/(1+1e-7*t*k))*(inv(W)'+ (eye(size(y))-tanh(y))*x');
        W=dW+W;
    end
    W = normRows(W);
%     delta = max(1 - abs(dot(W,Wlast,2))); %stop iterations if converged
end

% save files and correlation

unmixed = rescale(My*W',-1,1);
U1files = {'1unmixed1.wav' '1unmixed2.wav' '1unmixed3.wav'};
for i=1:length(U1files)
    audiowrite(U1files{i},unmixed(:,i),Mfs(i));
end
cormat1 = corr(unmixed, Sy)
save('table1','cormat1');
%% Q2
clear My
clear Mfs
clear Mfiles
Mfiles = {'noisy_mix1.wav' 'noisy_mix2.wav' 'noisy_mix3.wav'};
for i=1:length(Mfiles)
    [My(:,i),Mfs(:,i)] = audioread(Mfiles{i});
end

%initialize W
normRows = @(X) bsxfun(@rdivide,X,sqrt(sum(X.^2,2)));
W = normRows(rand(r,size(My,2))); % Random initial weights
delta=1;
k=0;

while delta > conLim && k < iterMax
    k=k+1;
    Wlast=W;
    sample = randperm(length(My));
    for t=1:length(My)
        x = My(sample(t),:)';
        y=W*x; 
        dW=(eta/(1+1e-7*t*k))*(inv(W)'+ (eye(size(y))-tanh(y))*x');
        W=dW+W;
    end
    W = normRows(W);
%     delta = max(1 - abs(dot(W,Wlast,2))); %stop iterations if converged
end

% save files and correlation

unmixedNoise = rescale(My*W',-1,1);
U2files = {'2unmixed1.wav' '2unmixed2.wav' '2unmixed3.wav'};
for i=1:length(U2files)
    audiowrite(U2files{i},unmixedNoise(:,i),Mfs(i));
end
cormat2 = corr(unmixedNoise, Sy)
save('table2','cormat2');

%% Q3
%sanger
clear My
clear Mfs
clear Mfiles
Mfiles = {'noisy_mix1.wav' 'noisy_mix2.wav' 'noisy_mix3.wav' ...
    'noisy_mix4.wav' 'noisy_mix5.wav' 'noisy_mix6.wav' 'noisy_mix7.wav'};
for i=1:length(Mfiles)
    [My(:,i),Mfs(:,i)] = audioread(Mfiles{i});
end
%initialize W
normRows = @(X) bsxfun(@rdivide,X,sqrt(sum(X.^2,2)));
W = normRows(rand(r,size(My,2))); % Random initial weights
delta=1;
k=0;
while delta > conLim && k < iterMax
    k=k+1;
    Wlast=W;
    sample = randperm(length(My));
    for t=1:length(My)
        x = My(sample(t),:)';
        y=W*x; 
        dW=(eta/(1+1e-7*t*k))*(y*x'-tril(y*y')*W);
        W=dW+W;
    end
    W = normRows(W);
%     delta = max(1 - abs(dot(W,Wlast,2))); %stop iterations if converged
end
save('eigenSanger','W');


%pca
coef=pca(My,'NumComponents',3 );
a=corrcoef(coef,W');
sangVSpca=a(1,2);
save('SangVSpca','sangVSpca')

%use infomax on sanger
Xsanger=W';
compressed_data=(Xsanger'*My')';
r=size(Sy,2); % #of independent components
eta=0.3; %learning rate
iterMax=10; %max iterations
conLim=1e-5; %convergence criteria
%initialize W
normRows = @(X) bsxfun(@rdivide,X,sqrt(sum(X.^2,2)));
W = normRows(rand(r,size(compressed_data,2))); % Random initial weights
delta=1;
k=0;
while delta > conLim && k < iterMax
    k=k+1;
    Wlast=W;
    sample = randperm(length(compressed_data));
    for t=1:length(compressed_data)
        x = compressed_data(sample(t),:)';
        y=W*x; 
        dW=(eta/(1+1e-7*t*k))*(inv(W)'+ (eye(size(y))-tanh(y))*x');
        W=dW+W;
    end
    W = normRows(W);
%     delta = max(1 - abs(dot(W,Wlast,2))); %stop iterations if converged
end
unmixedSanger = rescale(compressed_data*W',-1,1);
cormat3 = corr(unmixedSanger, Sy)
save('table3','cormat3');
U3files = {'3unmixed1.wav' '3unmixed2.wav' '3unmixed3.wav'};
for i=1:length(U3files)
    audiowrite(U3files{i},unmixedSanger(:,i),Mfs(i));
end