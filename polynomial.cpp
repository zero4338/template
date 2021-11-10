#include<bits/stdc++.h>
typedef long long ll;
using namespace std;
int read()
{
	int ret=0;bool f=0;char c=getchar();
	while(c>'9'||c<'0')f|=(c=='-'),c=getchar();
	while(c>='0'&&c<='9')ret=(ret<<3)+(ret<<1)+(c^48),c=getchar();
	return f?-ret:ret;
}
const int mod=998244353;
int qpow(int a,int b){int ret=1;for(;b;b>>=1,a=(ll)a*a%mod)if(b&1)ret=(ll)ret*a%mod;return ret;}
int n,k;
bool deb=0;
int R[1<<21],W[1<<21],iinv[1<<21];
struct poly
{
	vector<int>v;
	int& operator [](int i){return v[i];}
	void ntt(int L,int typ)
	{
		int n=pow(2,L);
		for(int i=0;i<n;i++)R[i]=(R[i>>1]>>1)|((i&1)<<(L-1));
		W[0]=1;W[1]=qpow(3,(mod-1)/n);if(typ==-1)W[1]=qpow(W[1],mod-2);
		for(int i=2;i<n;i++)W[i]=(ll)W[i-1]*W[1]%mod;
		if(v.size()<n)v.resize(n,0);
		for(int i=0;i<n;i++)if(R[i]>i)swap(v[i],v[R[i]]);
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
	poly operator *(poly &x)
	{
		poly ret;
		int ori=v.size(),xori=x.v.size();
		int L=ceil(log2(v.size()+x.v.size()-1)),n=pow(2,L);
		ntt(L,1);x.ntt(L,1);
		ret.v.resize(n);
		for(int i=0;i<n;i++)ret[i]=(ll)v[i]*x[i]%mod;
		ret.ntt(L,-1);ntt(L,-1);x.ntt(L,-1);
		while(v.size()>1&&v.back()==0)v.pop_back();
		while(x.v.size()>1&&x.v.back()==0)x.v.pop_back();
		while(ret.v.size()>1&&ret.v.back()==0)ret.v.pop_back();
		return ret;
	}
	void operator *=(poly &x)
	{
		int L=ceil(log2(v.size()+x.v.size())),n=pow(2,L);
		ntt(L,1);x.ntt(L,1);
		for(int i=0;i<n;i++)v[i]=(ll)v[i]*x[i]%mod;
		ntt(L,-1);x.ntt(L,-1);
		while(v.size()>1&&v.back()==0)v.pop_back();
		while(x.v.size()>1&&x.v.back()==0)x.v.pop_back();
	}
	void getinv(poly &ret,int deg=-1)
	{
		if(deg==-1)deg=v.size();
		if(deg==1){ret.v.resize(1);ret[0]=qpow(v[0],mod-2);return;}
		getinv(ret,(deg+1)>>1);
		int L=ceil(log2(deg*2)),n=pow(2,L);
		poly tmp;tmp.v.resize(deg+1,0);
		for(int i=0;i<=deg;i++)if(i<v.size())tmp.v[i]=v[i];
		tmp.ntt(L,1);ret.ntt(L,1);
		for(int i=0;i<n;i++)ret[i]=(ll)ret[i]*(2+mod-(ll)tmp.v[i]*ret[i]%mod)%mod;
		ret.ntt(L,-1);ret.v.resize(deg);
	}
	void getln(poly &ret)
	{
		poly rev;getinv(rev);ret.v.resize(v.size(),0);
		for(int i=0;i+1<v.size();i++)ret[i]=(ll)v[i+1]*(i+1)%mod;
		ret*=rev;ret.v.resize(v.size());
		iinv[1]=1;
		for(int i=2;i<ret.v.size();i++)iinv[i]=(ll)(mod-mod/i)*iinv[mod%i]%mod;
		for(int i=ret.v.size()-1;i>=1;i--)ret.v[i]=(ll)ret.v[i-1]*iinv[i]%mod;
		ret.v[0]=0;
	}
	void getexp(poly &ret,int deg=-1)
	{
		if(deg==-1)deg=v.size();
		if(deg==1){ret.v.resize(1);ret[0]=1;return;}
		getexp(ret,(deg+1)>>1);
		ret.v.resize(deg);
		poly ln;ret.getln(ln);
		int L=ceil(log2(deg*2)),n=pow(2,L);
		ln[0]=(v[0]+1-ln[0]+mod)%mod;
		for(int i=1;i<deg;i++)ln[i]=(v[i]-ln[i]+mod)%mod;
		ret*=ln;
		ret.v.resize(deg);
	}
	void getpow(poly &ret,int k)
	{
		poly ln;getln(ln);
		for(int &i:ln.v)i=(ll)i*k%mod;
		ln.getexp(ret);
	}
};
poly f,g;
int main()
{
	n=read();k=read();f.v.resize(n);generate_n(f.v.begin(),n,read);
	f.getpow(g,k);
	for(int &i:g.v)printf("%d ",i);printf("\n");
	return 0;
}