function [score, label] = run_12ECG_classifier(data,header_data,classes, model)
    th = load('optimalthreshold');
    x=th.x;
    x=x/10000;
    num_classes = length(classes);

    label = zeros([1,num_classes]);
    score = ones([1,num_classes]);
    output2=[0 0 0 0 0 0 0 0 0];
    % Use your classifier here to obtain a label and score for each class.
   [fb_v1, fb_avr, avb_feat, lbbb_feat, normal_feat, fb_ar] = get_12ECG_features(data,header_data);
    [labelaf, scoreaf]=predict(model{1,1},fb_avr);
    [labelavb, scoreavb]=predict(model{2,1},avb_feat);
    [labellbbb, scorelbbb]=predict(model{3,1},lbbb_feat);
    [labelnormal, scorenormal]=predict(model{4,1},normal_feat);
    [labelpac, scorepac]=predict(model{5,1},fb_ar);
    [labelpvc, scorepvc]=predict(model{6,1},fb_ar);
    [labelrbbb, scorerbbb]=predict(model{7,1},fb_v1);
    [labelstd, scorestd]=predict(model{8,1},fb_ar);
    [labelste, scoreste]=predict(model{9,1},fb_ar);
   score=[scoreaf(1,1) scoreavb(1,1) scorelbbb(1,1) scorenormal(1,1) scorepac(1,1) scorepvc(1,1) scorerbbb(1,1) scorestd(1,1) scoreste(1,1)];
output1=[scoreaf(1,1)<x(1) scoreavb(1,1)<x(2) scorelbbb(1,1)<x(3) scorenormal(1,1)<x(4) scorepac(1,1)<x(5) scorepvc(1,1)<x(6) scorerbbb(1,1)<x(7) scorestd(1,1)<x(8) scoreste(1,1)<x(9)];
output1=double(output1);
% if sum(output1)>3
% score(output1==0)=0;
% b=sort(score,'descend');
% max3=b(1:3);
% id(1)=find((score)==b(1));
% id(2)=find((score)==b(2));
% id(3)=find((score)==b(3));
% output2(id)=1;
% label=output2;
% else
%     label=output1;
% end
if sum(output1)==0
    a=find(max(score));
    output2(a)=1;
    label=output2;
else
    label=output1;
end
    
end



