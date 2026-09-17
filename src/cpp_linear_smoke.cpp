#include "openfhe.h"
#include <complex>
#include <iostream>
#include <vector>
#include <cmath>
using namespace lbcrypto;
int main(){
 using C=std::complex<double>; CCParams<CryptoContextCKKSRNS> p; p.SetMultiplicativeDepth(5); p.SetScalingModSize(50); p.SetRingDim(1<<15); p.SetSecurityLevel(HEStd_128_classic); p.SetScalingTechnique(FIXEDMANUAL); p.SetCKKSDataType(COMPLEX); auto cc=GenCryptoContext(p); for(auto f:{PKE,KEYSWITCH,LEVELEDSHE,ADVANCEDSHE})cc->Enable(f);
 auto kp=cc->KeyGen(); cc->EvalMultKeyGen(kp.secretKey); std::vector<int32_t> keys; for(int i=1;i<8;i++) keys.push_back(i); cc->EvalRotateKeyGen(kp.secretKey,keys);
 std::vector<std::vector<C>> A(8,std::vector<C>(8)); for(int i=0;i<8;i++)for(int j=0;j<8;j++)A[i][j]=(i==j?1.0:0.01*(i+j));
 auto pre=cc->EvalLinearTransformPrecompute(A); std::vector<C>x(cc->GetRingDimension()/2); for(int i=0;i<8;i++)x[i]=C(i+1,10+i); auto ct=cc->Encrypt(kp.publicKey,cc->MakeCKKSPackedPlaintext(x)); auto out=cc->EvalLinearTransform(pre,ct); Plaintext d; cc->Decrypt(out,kp.secretKey,&d); d->SetLength(8); auto y=d->GetCKKSPackedValue(); for(int i=0;i<8;i++){C e=0;for(int j=0;j<8;j++)e+=A[i][j]*x[j];std::cout<<i<<" "<<y[i]<<" exp "<<e<<"\n";} }
