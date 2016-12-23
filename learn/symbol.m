mupad;
r=10;
f=sym('x*exp(-x^10)*10*x^(10-1)');
ezplot(f,[-2 2]);
subexpr(x,10)
limit(f);
limit(f,x,1)
limit(f,x,inf);
int(f);%%Çó»ý·Ö



%% help whittakerM;
ezplot(sym('whittakerM(-1/(2*1), 1/(2*1) + 1/2, -x^1)'),[-100 100])