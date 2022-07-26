class poly
{
	static int R[1<<21];
	static mint W[1<<21];
	private:
	vector<mint>v;
	void ntt(int L,int typ)
	{
		int n=1<<L;
		for(int i=0;i<n;i++)R[i]=(R[i>>1]>>1)|((i&1)<<(L-1));
		W[0]=1;W[1]=qpow((mint)3,(mod-1)/n);if(typ==-1)W[1]=qpow(W[1],mod-2);
		for(int i=2;i<n;i++)W[i]=W[i-1]*W[1];
		set(n);
		for(int i=0;i<n;i++)if(R[i]>i)swap(v[R[i]],v[i]);
		for(int t=n>>1,d=1;d<n;d<<=1,t>>=1)
			for(int i=0;i<n;i+=(d<<1))
				for(int j=0;j<d;j++)
				{
					mint tmp=W[t*j]*v[i+j+d];
					v[i+j+d]=v[i+j]-tmp;
					v[i+j]+=tmp;
				}
		if(typ==-1){mint inv=qpow((mint)n,mod-2);for(int i=0;i<n;i++)v[i]*=inv;}
	}
	void adjust(){while(len()>1&&v.back().v==0)v.pop_back();}
	public:
	poly()=default;
	poly(initializer_list<mint>_v):v(_v){}
	void init(int n){v.resize(n);generate(v.begin(),v.end(),read);}
	mint&operator[](const int &i){return v[i];}
	const mint& operator[](const int &i)const{return v[i];}
	int len()const{return v.size();}
	void set(int l){v.resize(l);}
	poly operator *(poly rhs)
	{
		poly ret,tmp=*this;
		int L=ceil(log2(tmp.len()+rhs.len()-1)),n=1<<L;
		tmp.ntt(L,1);rhs.ntt(L,1);
		ret.set(n);
		for(int i=0;i<n;i++)ret[i]=tmp[i]*rhs[i];
		ret.ntt(L,-1);ret.adjust();
		return ret;
	}
	void operator *=(poly rhs)
	{
		int L=ceil(log2(len()+rhs.len()-1)),n=1<<L;
		ntt(L,1);rhs.ntt(L,1);
		for(int i=0;i<n;i++)v[i]*=rhs[i];
		ntt(L,-1);adjust();
	}
	poly operator +(const poly &rhs)
	{
		poly ret;
		ret.set(max(len(),rhs.len()));
		for(int i=0;i<len();i++)ret[i]=v[i];
		for(int i=0;i<rhs.len();i++)ret[i]+=rhs[i];
		return ret;
	}
	void operator +=(const poly &rhs)
	{
		if(len()<rhs.len())set(rhs.len());
		for(int i=0;i<rhs.len();i++)v[i]+=rhs[i];
	}
	poly operator -(const poly &rhs)
	{
		poly ret;
		ret.set(max(len(),rhs.len()));
		for(int i=0;i<len();i++)ret[i]=v[i];
		for(int i=0;i<rhs.len();i++)ret[i]-=rhs[i];
		return ret;
	}
	void operator -=(const poly &rhs)
	{
		if(len()<rhs.len())set(rhs.len());
		for(int i=0;i<rhs.len();i++)v[i]-=rhs[i];
	}	
	pair<poly,poly>operator %(poly rhs)
	{
		if(rhs.len()>len())return {{{0}},*this};
		poly tmp0=*this,tmp1=rhs;
		reverse(tmp0.v.begin(),tmp0.v.end());
		reverse(tmp1.v.begin(),tmp1.v.end());
		tmp1=tmp1.getinv(len()-rhs.len()+1);
		tmp0*=tmp1;tmp0.set(len()-rhs.len()+1);
		reverse(tmp0.v.begin(),tmp0.v.end());
		poly r=*this;
		tmp1=tmp0*rhs;r=r-tmp1;r.set(rhs.len()-1);
		return {tmp0,r};
	}
	poly getinv(int deg=-1)
	{
		if(deg==-1)deg=len();
		if(deg==1)return {{qpow(v[0],mod-2)}};
		poly ret=getinv((deg+1)>>1);
		int L=ceil(log2(deg))+1,n=1<<L;
		poly tmp;tmp.set(deg);
		for(int i=0;i<min(len(),deg);i++)tmp[i]=v[i];
		tmp.ntt(L,1);ret.ntt(L,1);
		for(int i=0;i<n;i++)ret[i]=ret[i]*(2-tmp[i]*ret[i]);
		ret.ntt(L,-1);ret.v.resize(deg);
		return ret;
	}
	poly getln(int deg=-1)
	{
		if(deg==-1)deg=len();
		poly ret,rev=getinv();
		ret.set(deg);
		for(int i=0;i+1<len();i++)ret[i]=v[i+1]*(i+1);
		ret*=rev;
		ret.set(deg);
		static mint inv[1<<21];
		inv[1]=1;for(int i=2;i<ret.len();i++)inv[i]=(mod-mod/i)*inv[mod%i];
		for(int i=ret.len()-1;i>=1;i--)ret[i]=ret[i-1]*inv[i];
		ret[0]=0;
		return ret;
	}
	poly getexp(int deg=-1)
	{
		if(deg==-1)deg=len();
		if(deg==1)return {{1}};
		poly ret=getexp((deg+1)>>1);
		ret.set(deg);
		poly ln=ret.getln();
		int L=ceil(log2(deg))+1,n=1<<L;
		ln[0]=v[0]+1-ln[0];
		for(int i=1;i<deg;i++)ln[i]=v[i]-ln[i];
		ret*=ln;ret.set(deg);
		return ret;
	}
	poly getpow(int k)
	{
		poly ln=getln();
		for(mint &i:ln.v)i*=k;
		return ln.getexp();
	}
	vector<mint>meval(vector<mint>x)
	{
		vector<poly>tmp(x.size()<<2);
		function<void(int,int,int)>get=[&](int u,int l,int r)
		{
			if(l==r){tmp[u]=poly{{-x[l],1}};return;}
			int mid=(l+r)>>1;
			get(u<<1,l,mid);get(u<<1|1,mid+1,r);
			tmp[u]=tmp[u<<1]*tmp[u<<1|1];
		};
		get(1,0,x.size()-1);
		vector<mint>ret;ret.resize(x.size());
		function<void(int,int,poly,int)>solve=[&](int l,int r,poly f,int u)
		{
			if(l==r){ret[l]=f[0];return;}
			int mid=(l+r)>>1;
			solve(l,mid,(f%tmp[u<<1]).second,u<<1);
			solve(mid+1,r,(f%tmp[u<<1|1]).second,u<<1|1);
		};
		solve(0,x.size()-1,(*this%tmp[1]).second,1);
		return ret;
	}
	void print(){for(mint &i:v)printf("%d ",i.v);putchar('\n');}
};
int poly::R[];
mint poly::W[];