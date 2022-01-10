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
int mod;
namespace Exlucas
{
	const int maxn=2e4+5;
	const int maxm=10;
	int mul[maxm],fac[maxm][maxn];
	vector<pair<int,int>>factor;
	int qpow(int a,ll b,int mod){int ret=1;for(;b;b>>=1,a=(ll)a*a%mod)if(b&1)ret=(ll)ret*a%mod;return ret;}
	void init(int mod)
	{
		int v=mod;
		for(int i=2;i<=v/i;i++)
			if(v%i==0)
			{
				factor.push_back({i,0});
				while(v%i==0)v/=i,factor.back().second++;
			}
		if(v!=1)factor.push_back({v,1});
		for(int i=0;i<factor.size();i++)mul[i]=qpow(factor[i].first,factor[i].second,mod+10);
		for(int i=0;i<factor.size();i++)
		{
			fac[i][0]=1;
			for(int j=1;j<=mul[i];j++)
			{
				fac[i][j]=fac[i][j-1];
				if(j%factor[i].first)fac[i][j]=(ll)fac[i][j]*j%mul[i];
			}
		}
		return;
	}
	void exgcd(int a,int b,int &x,int &y)
	{
		if(b==0){x=1;y=0;return;}
		exgcd(b,a%b,x,y);
		int z=x;x=y;y=z-y*(a/b);
	}
	int rev(int a,int b){int x0,y0;exgcd(a,b,x0,y0);x0=(x0%b+b)%b;return x0;}
	int exlucas(int n,int id)
	{
		if(n<=0)return 1;
		int ret=1;
		if(n>=mul[id])ret=qpow(fac[id][mul[id]],n/mul[id],mul[id]);
		for(int i=2;i<=n%mul[id];i++)if(i%factor[id].first)ret=(ll)ret*i%mul[id];
		ret=(ll)ret*exlucas(n/factor[id].first,id)%mul[id];
		return ret;
	}
	int C(int n,int m,int id)
	{
		int ret=1;
		int last=0;
		for(int i=n;i;i/=factor[id].first)last+=i/factor[id].first;
		for(int i=m;i;i/=factor[id].first)last-=i/factor[id].first;
		for(int i=n-m;i;i/=factor[id].first)last-=i/factor[id].first;
		if(last>=factor[id].second)return 0;
		return (ll)qpow(factor[id].first,last,mul[id])*exlucas(n,id)%mul[id]*rev(exlucas(m,id),mul[id])%mul[id]*rev(exlucas(n-m,id),mul[id])%mul[id];
	}
	int a[maxm];
	int getans()
	{
		static int M[maxm];
		int ret=0;
		for(int i=0;i<factor.size();i++)M[i]=mod/mul[i];
		for(int i=0;i<factor.size();i++)
		{
			int x0,y0;
			exgcd(M[i],mul[i],x0,y0);
			x0=(x0%mul[i]+mul[i])%mul[i];
			(ret+=(ll)M[i]*a[i]%mod*x0%mod)%=mod;
		}
		return ret;
	}
	int C(int n,int m)
	{
		if(n<m||m<0)return 0;
		for(int i=0;i<factor.size();i++)a[i]=C(n,m,i);
		return getans();
	}
}
using Exlucas::init;
using Exlucas::C;
int main()
{
	mod=read();init(mod);
	return 0;
}
