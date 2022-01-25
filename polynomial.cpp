#include<bits/stdc++.h>
typedef long long ll;
using namespace std;
const int mod=998244353;
int inv[1<<21];
void prework(){inv[1]=1;for(int i=2;i<(1<<21);i++)inv[i]=(ll)(mod-mod/i)*inv[mod%i]%mod;}
int qpow(int a,int b){int ret=1;for(;b;b>>=1,a=(ll)a*a%mod)if(b&1)ret=(ll)ret*a%mod;return ret;}
int R[1<<21],W[1<<21];
struct poly
{
	vector<int>v;
	int&operator[](const int &i){return v[i];}
	int len(){return v.size();}
	void set(int l){v.resize(l);}
	void ntt(int L,int typ)
	{
		int n=1<<L;
		for(int i=0;i<n;i++)R[i]=(R[i>>1]>>1)|((i&1)<<(L-1));
		W[0]=1;W[1]=qpow(3,(mod-1)/n);if(typ==-1)W[1]=qpow(W[1],mod-2);
		for(int i=2;i<n;i++)W[i]=(ll)W[i-1]*W[1]%mod;
		set(n);
		for(int i=0;i<n;i++)if(R[i]>i)swap(v[R[i]],v[i]);
		for(int t=n>>1,d=1;d<n;d<<=1,t>>=1)
			for(int i=0;i<n;i+=(d<<1))
				for(int j=0;j<d;j++)
				{
					int tmp=(ll)W[t*j]*v[i+j+d]%mod;
					v[i+j+d]=(v[i+j]-tmp+mod)%mod;
					v[i+j]=(v[i+j]+tmp)%mod;
				}
		if(typ==-1){int inv=qpow(n,mod-2);for(int i=0;i<n;i++)v[i]=(ll)v[i]*inv%mod;}
	}
	void adjust(){while(len()>1&&v.back()==0)v.pop_back();}
	poly operator *(poly &x)
	{
		poly ret,tmp0=*this,tmp1=x;
		int L=ceil(log2(tmp0.len()+tmp1.len()-1)),n=1<<L;
		tmp0.ntt(L,1);tmp1.ntt(L,1);
		ret.set(n);
		for(int i=0;i<n;i++)ret[i]=(ll)tmp0[i]*tmp1[i]%mod;
		ret.ntt(L,-1);ret.adjust();
		return ret;
	}
	poly operator +(poly &x)
	{
		poly ret;
		ret.set(max(len(),x.len()));
		for(int i=0;i<len();i++)ret[i]=v[i];
		for(int i=0;i<x.len();i++)ret[i]=(ret[i]+x[i])%mod;
		return ret;
	}
	void operator +=(poly &x)
	{
		if(len()<x.len())set(x.len());
		for(int i=0;i<x.len();i++)v[i]=(v[i]+x[i]+mod)%mod;
	}
	poly operator -(poly &x)
	{
		poly ret;
		ret.set(max(len(),x.len()));
		for(int i=0;i<len();i++)ret[i]=v[i];
		for(int i=0;i<x.len();i++)ret[i]=(ret[i]-x[i]+mod)%mod;
		return ret;
	}
	void operator -=(poly &x)
	{
		if(len()<x.len())set(x.len());
		for(int i=0;i<x.len();i++)v[i]=(v[i]-x[i]+mod)%mod;
	}
	void operator *=(poly &x)
	{
		poly tmp=x;
		int L=ceil(log2(len()+tmp.len()-1)),n=1<<L;
		ntt(L,1);tmp.ntt(L,1);
		for(int i=0;i<n;i++)v[i]=(ll)v[i]*tmp[i]%mod;
		ntt(L,-1);adjust();
	}
	pair<poly,poly>operator %(poly &x)
	{
		if(x.len()>len())return {{{0}},*this};
		poly tmp0=*this,tmp1=x;
		reverse(tmp0.v.begin(),tmp0.v.end());
		reverse(tmp1.v.begin(),tmp1.v.end());
		tmp1=tmp1.getinv(len()-x.len()+1);
		tmp0*=tmp1;tmp0.set(len()-x.len()+1);
		reverse(tmp0.v.begin(),tmp0.v.end());
		poly r=*this;
		tmp1=tmp0*x;r=r-tmp1;r.set(x.len()-1);
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
		for(int i=0;i<n;i++)ret[i]=(ll)ret[i]*(2-(ll)tmp[i]*ret[i]%mod+mod)%mod;
		ret.ntt(L,-1);ret.v.resize(deg);
		return ret;
	}
	poly getln(int deg=-1)// call prework before this
	{
		if(deg==-1)deg=len();
		poly ret,rev=getinv();
		ret.set(deg);
		for(int i=0;i+1<len();i++)ret[i]=(ll)v[i+1]*(i+1)%mod;
		ret=ret*rev;
		ret.set(deg);
		for(int i=ret.len()-1;i>=1;i--)ret[i]=(ll)ret[i-1]*inv[i]%mod;
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
		ln[0]=(v[0]+1-ln[0]+mod)%mod;
		for(int i=1;i<deg;i++)ln[i]=(v[i]-ln[i]+mod)%mod;
		ret=ret*ln;ret.set(deg);
		return ret;
	}
	poly getpow(int k)
	{
		poly ln=getln();
		for(int &i:ln.v)i=(ll)i*k%mod;
		return ln.getexp();
	}
	vector<int>multipointeval(vector<int>x)
	{
		static poly tmp[1<<21];
		auto get=[&](int u,int l,int r)
		{
			function<void(int,int,int)>s;
			s=[&](int u,int l,int r)
			{
				if(l==r){tmp[u]=poly{{mod-x[l],1}};return;}
				int mid=(l+r)>>1;
				s(u<<1,l,mid);s(u<<1|1,mid+1,r);
				tmp[u]=tmp[u<<1]*tmp[u<<1|1];	
			};
			return s(u,l,r);
		};
		get(1,0,x.size()-1);
		vector<int>ret;ret.resize(x.size());
		auto solve=[&](int l,int r,poly f,int u)
		{
			function<void(int,int,poly,int)>s;
			s=[&](int l,int r,poly f,int u)
			{
				if(l==r){ret[l]=f[0];return;}
				int mid=(l+r)>>1;
				s(l,mid,(f%tmp[u<<1]).second,u<<1);
				s(mid+1,r,(f%tmp[u<<1|1]).second,u<<1|1);
			};
			return s(l,r,f,u);
		};
		solve(0,x.size()-1,(*this%tmp[1]).second,1);
		return ret;
	}
	void print(){for(int &i:v)printf("%d ",i);putchar('\n');}
};
int main()
{
	
}