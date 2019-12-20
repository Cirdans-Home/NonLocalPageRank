function xx = linearize(x,a,b)
%% linearize.m
% Code by:
% S. Cipolla - Università di Padova, Dipartimento di Matematica
% F. Durastante - Consiglio Nazionale delle Ricerche, Istituto per le
% Applicazioni del Calcolo "M. Picone"
% F. Tudisco - Gran Sasso Science Institute
    xx = x;
    alpha = (b-a)./(max(max(x))-min(min(x(x>0))));
    beta = alpha*min(min(x(x>0)))-a;
    xx(x>0) = alpha.*x(x>0)-beta;
end