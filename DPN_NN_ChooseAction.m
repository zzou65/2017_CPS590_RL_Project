
function actionidx = DPN_NN_ChooseAction(nn, naction, if_max_prob_act)
        feature = AD2D_Generate_Feature();
        feature = feature(:)';
        feature = normalize(feature, nn.statemu, nn.statesigma);
        nn = nnff(nn, feature, zeros(1, nn.size(end)));
        nnoutvec = nn.a{end}(:);
        expout = exp(nnoutvec(1:naction));
        actionprobs = expout./sum(expout);
        if(if_max_prob_act == 1)
            [~, actionidx] = max(actionprobs);
            return;
        end
        probsinterval = cumsum(actionprobs);
        probsinterval = [0; probsinterval(:)];
        if(abs(probsinterval(end)-1) > 1e-10)
            error('Check me: total action prob not 1');
        end
        rd = rand(1);
        actionidx = 0;
        for i=1:naction
            if(rd > probsinterval(i) && rd < probsinterval(i+1))
                actionidx = i;
                break;
            end
        end
        if(actionidx == 0)
            [~, actionidx] = max(actionprobs);
        end
end
    


