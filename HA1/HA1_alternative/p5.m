function val = p5(x1,x2,y,sigma,pos_vec)

    v=90; gamma=3; %Constatns
    xs=zeros(12,length(x1)); % Placeholders
    xs(1:6,:)=bsxfun(@minus,x1,transpose(pos_vec(1,:)));
    xs(7:12,:)=bsxfun(@minus,x2,transpose(pos_vec(2,:)));


    val=prod(normpdf(y,v-10*gamma*log10(sqrt(xs(1:6,:).^2+xs(7:12,:).^2)),sigma));
end

