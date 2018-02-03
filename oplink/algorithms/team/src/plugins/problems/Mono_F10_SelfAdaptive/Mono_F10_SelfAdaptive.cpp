/***********************************************************************************
 * AUTORES
 *   - Eduardo Manuel Segredo Gonz�lez
 *   - Carlos Segura Gonz�lez
 * 
 * FECHA
 *    Abril 2012
 ************************************************************************************/

#include <math.h>
#include "Mono_F10_SelfAdaptive.h"

double f10[1000] = {-6.571740e+00,-2.569351e+00,-4.693453e+00,-1.221222e+00,4.608476e+00,-4.172883e+00,4.299095e+00,-5.564266e+00,-3.154860e+00,5.014815e+00,-1.948804e+00,-6.029759e+00,-6.757910e+00,-4.969640e+00,6.335886e+00,1.799810e+00,6.408511e+00,3.985204e+00,1.909531e+00,-1.559749e-02,4.368818e+00,2.982258e+00,-2.024213e+00,6.049229e+00,3.383029e+00,6.557614e-02,7.185533e-01,5.575540e+00,-5.122889e+00,-9.557999e-01,5.136575e+00,-6.108188e+00,-7.941781e-01,-1.479978e+00,-1.313630e+00,-6.460725e+00,2.247470e-01,-5.276668e+00,-6.902121e+00,3.627432e+00,-6.541585e+00,-1.165224e+00,-4.688885e+00,-6.230674e+00,3.349402e+00,-3.306404e+00,-4.103180e+00,2.521686e+00,2.249735e+00,-6.072993e+00,
-1.281808e+00,4.593614e+00,2.941809e+00,3.486103e+00,1.535505e+00,3.451281e+00,-1.198811e+00,3.616844e+00,3.015630e+00,6.441931e-02,1.410862e+00,2.793187e+00,-2.290243e+00,-2.090054e+00,-4.471412e-02,-3.837954e+00,-3.418125e+00,1.752517e+00,-6.253379e+00,-5.058774e+00,-6.595191e+00,-7.802374e-01,4.309114e+00,-2.497726e+00,-3.298718e+00,-1.473184e-01,1.912453e+00,6.316808e+00,-6.157324e+00,-4.492617e+00,-3.142174e+00,6.089071e+00,4.879015e+00,-3.360108e+00,2.796324e+00,2.552950e-01,2.626477e-01,-3.295541e+00,3.862188e+00,2.095133e+00,-6.960492e+00,-3.735683e+00,2.799186e+00,-6.292280e+00,-5.037049e+00,5.043402e+00,-2.370358e+00,-5.308549e+00,4.735217e+00,-6.200087e+00,
3.796112e+00,6.153834e+00,-2.824421e+00,-3.319247e+00,1.891198e+00,-5.816845e+00,-4.252645e+00,-3.530729e+00,-5.152457e+00,3.661117e+00,1.182498e+00,-4.134473e+00,2.984709e+00,1.250515e+00,4.268547e+00,-1.843671e+00,5.307079e+00,-3.826425e+00,2.315439e+00,-5.731888e+00,4.309682e+00,-1.804790e+00,-5.384292e+00,-1.232908e+00,-4.669317e+00,-5.432663e+00,-1.580425e+00,6.961001e+00,3.773289e+00,-2.460081e+00,-2.495726e+00,3.746192e+00,6.534252e+00,4.498323e+00,-2.828793e+00,-6.748989e+00,3.111254e+00,2.882948e+00,-2.611389e+00,-5.845766e+00,4.652847e+00,1.671618e+00,-2.975484e-01,1.920903e+00,-3.581002e+00,-2.353006e+00,-6.967197e+00,-2.783403e+00,-2.132249e+00,-3.055628e+00,
1.952579e-01,-2.366290e-01,-9.695837e-01,2.928084e+00,-1.564563e+00,6.498900e+00,6.608778e+00,-5.332431e+00,3.035570e+00,1.399549e+00,3.313541e+00,1.673422e+00,1.006641e+00,-4.273860e+00,1.748419e+00,-8.064325e-01,-4.991306e+00,-3.337187e+00,2.839283e+00,-3.615854e-01,2.074643e+00,2.652329e+00,3.096450e+00,-3.551768e+00,-3.243563e+00,-4.239754e-01,1.668520e+00,3.470957e-01,2.248853e+00,5.714212e-01,6.617582e+00,-3.689077e+00,1.830632e+00,1.507173e+00,2.243696e+00,-2.936985e+00,4.991149e+00,1.360962e+00,3.742487e+00,-6.612994e+00,-2.754580e+00,-5.357763e+00,-3.320286e+00,-6.387041e+00,-3.996556e+00,4.721609e+00,5.806434e+00,-1.908217e+00,-2.106093e+00,-7.462778e-01,
-7.978054e-01,-5.908783e+00,5.809159e+00,3.772623e+00,-6.309612e+00,4.346613e-01,-1.970235e+00,-3.520004e+00,-9.287028e-01,6.336180e+00,2.147758e+00,2.440160e+00,5.676654e+00,-3.629981e+00,5.605951e+00,3.205368e+00,-6.731715e+00,-4.132277e+00,-2.537744e+00,7.656301e-01,5.196495e+00,5.720908e+00,-3.439301e+00,-6.240772e+00,3.547278e+00,6.757263e+00,5.857844e+00,-6.666600e+00,1.027954e+00,3.716439e-01,4.649867e+00,-6.729833e+00,1.561808e-01,-6.247203e+00,-5.260826e+00,1.350512e+00,2.734796e+00,4.601790e+00,2.488668e+00,6.277633e+00,-4.511499e+00,1.447743e+00,1.878071e-01,5.821061e+00,2.637368e+00,1.987283e-01,3.838003e-01,-5.006933e+00,2.560508e+00,3.723840e+00,
1.873140e+00,3.579963e+00,5.670733e+00,1.224516e+00,3.803102e-01,5.242669e+00,-3.633471e+00,6.491410e+00,4.868486e+00,-2.074741e+00,5.676282e+00,2.801765e-01,-5.031402e+00,3.444870e+00,-3.239563e+00,-4.312741e+00,-3.759966e-01,-1.245907e+00,-1.452626e+00,-4.164982e+00,-6.985148e-01,2.785050e+00,-4.741099e+00,4.642279e+00,1.617679e+00,-3.617275e+00,-1.475429e+00,8.268044e-01,-3.853757e+00,1.561210e+00,-5.087714e+00,-4.010644e-02,3.632765e+00,-2.210599e+00,3.438380e+00,-2.826891e+00,-2.847028e+00,-2.892742e-01,2.279440e+00,5.030706e-01,2.301390e-01,-1.995607e+00,-2.488491e+00,3.257268e+00,-5.152829e+00,4.764196e+00,2.045252e+00,-6.843162e+00,-6.521330e+00,-6.913189e-01,
-1.242800e-01,-6.263379e+00,-4.472324e+00,-5.025971e+00,-4.612064e+00,-3.280767e-01,2.374025e+00,4.876270e+00,6.181225e+00,-3.536592e+00,4.141159e+00,3.780946e-01,6.943924e+00,-2.544234e+00,-6.212655e+00,2.165228e+00,4.773294e+00,-5.267727e+00,1.485193e+00,4.482010e+00,-6.450196e+00,-4.358896e+00,6.507909e-01,-6.528252e+00,2.358966e+00,4.761118e+00,-5.493621e+00,-4.194334e+00,5.373292e+00,5.015815e+00,-1.081423e+00,1.673265e+00,2.049056e+00,-3.757016e+00,6.206165e+00,5.242532e+00,2.328379e+00,5.163398e+00,6.947727e+00,1.155930e+00,6.795546e-01,6.825634e+00,1.795329e-01,7.242197e-01,1.355178e+00,-6.149853e+00,-3.991968e+00,1.705774e+00,4.543478e+00,-2.721522e+00,
1.439459e-01,-1.299072e-01,6.389610e+00,3.911285e+00,3.335638e+00,-5.202396e+00,-4.328290e+00,-1.693362e+00,1.365491e+00,1.259270e-01,9.010763e-01,-4.385346e+00,4.752785e+00,-6.840319e+00,1.875317e+00,-5.881529e+00,-8.194713e-01,-2.935613e+00,6.809967e+00,-8.682931e-01,6.917709e+00,-4.282723e+00,-3.146547e+00,-1.574896e+00,6.032877e+00,4.808175e+00,-4.224274e+00,4.046819e-01,2.125328e+00,2.930123e+00,-6.739597e+00,-6.399315e+00,6.668149e+00,5.746887e+00,2.229696e+00,1.138882e-01,5.823473e+00,-5.803650e+00,6.519821e+00,-3.477143e+00,6.037926e-01,2.036448e+00,3.496838e-01,3.358049e+00,6.063268e+00,6.370855e-01,-5.138742e-01,2.207991e+00,6.172970e+00,1.924589e+00,
4.494009e+00,3.623667e+00,-5.224141e+00,4.382719e+00,-1.966372e+00,5.053558e+00,5.911313e+00,3.200486e+00,2.647290e+00,-2.349653e+00,3.251582e+00,6.181029e+00,6.272594e+00,-5.136673e+00,-6.977413e+00,-3.801730e-01,2.938927e+00,-2.470257e+00,-2.880791e+00,-4.451325e+00,-3.398999e-01,3.040119e+00,1.421411e+00,5.514876e+00,5.062166e+00,-3.584649e+00,3.257386e+00,2.407690e+00,-6.292456e+00,3.069510e+00,2.298165e+00,-5.445976e+00,4.082318e+00,-9.455649e-01,6.811556e+00,5.449192e+00,-4.644181e+00,-4.914289e+00,3.931205e+00,-1.140480e+00,-5.284609e+00,-3.411341e+00,-3.338560e+00,-1.933766e+00,-2.638741e+00,6.959041e+00,-5.047637e+00,-6.344709e+00,3.406420e+00,6.756773e+00,
-1.312209e-01,-5.171280e+00,5.311059e+00,2.753992e+00,-5.807924e+00,9.745639e-01,5.287727e+00,6.969648e+00,9.188894e-02,2.357074e-01,-1.304376e+00,6.532114e+00,-2.787187e+00,3.887952e+00,4.806665e+00,-1.782653e+00,6.262222e+00,-6.736931e+00,-4.412071e+00,3.028766e+00,1.128147e+00,-2.240559e+00,-3.654137e+00,7.579440e-01,5.052470e-01,-3.181526e+00,4.752667e+00,2.565645e+00,1.281396e+00,-3.846365e+00,-2.424748e+00,-5.470495e-01,3.994184e+00,3.283405e+00,-4.689924e+00,-1.142205e+00,-4.137061e+00,-3.145860e+00,2.686798e+00,-3.502162e+00,3.429105e+00,-5.276021e+00,7.331410e-01,-9.630055e-02,-4.041398e+00,3.448124e+00,5.400105e-01,4.411953e+00,-2.134112e+00,-2.175708e-01,
3.253543e+00,-5.635518e+00,-4.276625e+00,2.886379e+00,-5.306373e+00,5.612852e+00,-1.320611e+00,-1.328453e+00,-6.149775e+00,-2.526382e-02,-4.120924e+00,-2.310301e+00,-2.629604e+00,-3.102391e+00,-1.456959e+00,-3.218407e+00,4.086710e+00,-5.238356e+00,5.409232e+00,-4.223372e+00,2.442003e+00,1.522290e+00,1.752145e+00,2.096387e+00,6.013407e+00,-2.756325e+00,-4.476716e+00,-2.452169e-01,6.029544e+00,1.446861e+00,-7.486110e-01,5.706026e+00,-4.065142e+00,-4.391562e+00,5.061489e-01,-4.038868e+00,5.095596e+00,-5.992094e+00,6.200126e+00,2.816382e+00,1.357962e+00,4.844605e+00,3.645157e+00,6.762126e+00,-1.086384e+00,-2.947936e-02,-5.888372e+00,2.775226e+00,6.404158e+00,1.990411e+00,
-5.993604e+00,-2.866664e-01,-3.120685e+00,3.527651e+00,3.084921e+00,-3.007787e+00,-6.902239e+00,-1.134715e+00,-3.641588e+00,1.919433e+00,-3.677959e+00,2.556352e+00,-3.724095e+00,-2.373378e+00,-2.657652e-01,4.927386e+00,-7.678064e-01,1.240672e+00,8.120598e-01,-1.177969e+00,6.553878e+00,1.291043e+00,-4.847879e+00,3.706331e+00,6.332474e+00,-4.113807e+00,1.195988e+00,3.321678e+00,3.153125e-01,4.221911e-01,-3.286856e+00,4.116513e+00,-6.452294e+00,-4.226058e+00,2.638270e+00,-5.747926e+00,5.637058e-03,5.549051e+00,-1.289239e+00,-4.633965e+00,-2.159581e+00,-5.424545e+00,-4.483657e+00,-3.565944e+00,-6.991000e+00,-2.382171e-01,-3.963538e+00,5.688771e+00,1.095089e+00,1.835828e+00,
1.934883e+00,1.609601e+00,-4.606986e+00,3.176850e-01,3.317698e+00,-2.300616e+00,-9.218305e-02,-2.941691e+00,4.688924e+00,4.467755e+00,-5.135712e+00,-3.189114e+00,6.768185e+00,2.224304e+00,6.573564e+00,2.909398e+00,-3.675803e+00,5.287433e+00,-6.370983e+00,-6.006760e+00,-5.466348e+00,3.423223e+00,1.298651e+00,-3.882501e+00,2.463355e+00,5.677478e+00,9.013312e-01,-4.244587e+00,2.606350e+00,2.324919e-01,4.686767e+00,-2.549842e+00,-3.076765e+00,-1.831769e+00,5.114890e+00,-6.741029e+00,-6.863985e+00,-5.408683e+00,2.934946e+00,9.429376e-01,3.658725e+00,-2.865106e+00,-4.159962e+00,-4.354848e-01,-6.241938e-02,6.855554e+00,6.917493e+00,-5.167280e+00,-1.875983e+00,1.999018e+00,
4.498225e+00,-1.125794e+00,-6.146893e+00,-6.037955e+00,1.638129e+00,9.937005e-01,-3.921872e+00,5.900548e+00,-6.720255e-01,-1.660912e+00,6.302034e-01,-6.622170e+00,-1.776506e-01,-1.814142e+00,-4.570164e+00,-4.646612e+00,-6.312475e+00,-5.960703e+00,1.961902e+00,-2.884850e+00,1.783575e+00,2.815970e+00,-6.805183e+00,-4.642553e+00,-2.882419e+00,1.687147e+00,4.601045e+00,-2.937123e+00,-7.638458e-01,5.318020e+00,-6.704187e+00,6.056072e+00,3.102999e+00,5.915234e+00,-6.870161e+00,-3.769074e+00,-1.708675e+00,-5.626617e+00,6.675364e+00,-4.997815e+00,-4.521204e+00,-3.918500e+00,-5.741456e+00,5.542384e+00,-5.354626e+00,-5.964595e-01,6.404688e+00,3.171801e+00,1.971686e+00,-3.152468e+00,
-4.634730e+00,4.522165e+00,3.253445e+00,8.259515e-02,-6.545781e+00,-1.501722e+00,3.500417e+00,5.756867e+00,-3.283209e+00,-5.057146e+00,-5.873765e+00,2.839597e+00,-5.014324e+00,4.337740e+00,-1.187047e+00,2.164248e+00,-4.745942e+00,3.341942e-01,-2.110465e+00,-4.083847e+00,-5.925842e+00,-6.882798e-01,6.468077e+00,-6.968089e-01,-3.059981e+00,1.644246e+00,-2.965337e+00,6.001604e+00,-1.058777e+00,1.089246e+00,-4.823096e+00,-4.680905e+00,3.257857e+00,-1.808476e+00,3.657195e+00,6.806919e-01,5.115007e+00,-3.382578e+00,3.296032e+00,-5.828041e+00,-1.316140e+00,3.539062e+00,1.137852e+00,-2.388397e+00,-3.179761e+00,-2.772109e+00,-1.629227e+00,6.101522e+00,5.126291e-02,6.212479e+00,
-3.627844e+00,-3.988772e+00,5.597794e+00,-5.220768e+00,-1.968176e+00,-1.547250e+00,-2.826274e-01,-3.027051e-01,4.201823e+00,3.896442e+00,-1.235280e+00,4.812410e+00,-1.875317e+00,8.544505e-01,4.389170e+00,9.947691e-02,3.580433e+00,3.710468e+00,-5.344646e+00,-2.411239e+00,1.528250e+00,4.149002e+00,2.981827e+00,5.669596e+00,-2.052585e+00,-6.008858e+00,-1.649285e+00,5.858158e+00,-1.373903e+00,-3.430733e+00,-5.961507e+00,-1.454626e+00,2.615271e+00,3.086333e+00,1.104755e+00,-3.542043e+00,6.129815e+00,1.889532e+00,-3.327207e+00,7.887861e-01,-3.735232e+00,9.631329e-01,5.307687e+00,-2.647701e+00,1.848965e+00,-5.235023e+00,6.127403e+00,-2.532353e+00,1.562396e-01,-4.995580e+00,
-4.588477e+00,-3.998380e+00,-2.714503e+00,-6.700030e+00,5.198377e+00,4.755314e+00,2.407180e+00,6.107482e+00,-2.725189e+00,-6.967609e+00,-4.792822e+00,-4.591476e+00,-2.253549e-01,-4.350642e+00,-1.457233e+00,3.200966e-01,-1.573572e-01,-5.446780e+00,-6.357395e+00,6.574584e+00,-6.909101e+00,-4.409012e+00,-5.604382e+00,2.970082e+00,-4.487519e+00,-6.122619e+00,-1.734616e+00,-4.401895e+00,1.338159e+00,5.838903e+00,6.790145e+00,-1.795790e+00,-5.444809e-01,6.205891e+00,6.223802e-01,-3.999929e+00,-3.876227e+00,-1.627580e+00,6.473411e+00,2.024331e+00,6.088777e+00,-5.002580e+00,-1.565288e+00,-6.734480e+00,6.391335e+00,3.204996e+00,5.982526e+00,4.637377e+00,6.147922e-01,3.874031e+00,
-4.865859e+00,-1.538289e+00,-2.290616e+00,-3.216525e+00,-2.451434e+00,3.429243e+00,-4.225244e-01,-2.616212e+00,-2.649525e+00,4.585810e+00,-3.575473e+00,3.629216e+00,3.996242e+00,-1.270240e+00,-8.049424e-01,5.554168e+00,4.111003e+00,-1.624855e+00,-4.158570e+00,-3.678312e+00,-4.711306e-01,3.911147e+00,2.262519e+00,-2.993003e+00,-6.505057e+00,8.306082e-01,-5.975310e+00,6.214028e+00,3.890128e+00,3.440282e+00,-1.467596e-02,1.059051e+00,6.779968e+00,-1.402863e+00,3.883737e+00,-5.365793e-01,3.305953e+00,-4.559674e+00,-1.660324e+00,-2.351869e+00,-6.482224e-01,-5.678301e+00,4.461010e+00,-5.809493e+00,-2.555744e+00,2.174992e+00,-1.189341e+00,5.843962e+00,5.987290e+00,3.167252e+00,
-2.536450e+00,1.362207e-01,2.707523e+00,1.916492e+00,-6.181646e-01,6.711853e+00,-4.559458e+00,-2.881183e+00,-4.235352e+00,-1.184890e+00,3.146174e+00,2.164973e+00,3.600854e-02,-3.998360e+00,-4.115121e+00,1.185321e+00,-1.992920e+00,3.091558e-01,-4.857350e+00,6.559212e+00,-5.341019e+00,-2.187796e+00,5.008521e+00,3.631647e+00,-3.605844e+00,1.398922e+00,6.655404e+00,5.702202e+00,-4.367602e+00,-5.802454e+00,5.434976e+00,1.130676e+00,-6.664639e+00,-1.795359e+00,-5.738064e+00,5.668919e-01,-2.137406e+00,-1.111520e+00,5.107321e+00,3.622314e+00,2.039899e+00,-5.873275e+00,-5.191867e+00,5.545620e+00,3.328609e-01,1.784732e+00,2.545136e+00,-7.894920e-01,3.194761e+00,5.592255e-01,
-1.308846e+00,2.406749e+00,4.652298e+00,-1.001779e+00,-4.016359e+00,-1.941138e+00,-5.833041e+00,-6.854515e+00,2.977219e+00,-3.011953e-01,-4.044378e+00,2.847204e+00,-6.231733e+00,-1.679069e+00,-2.728267e+00,-5.565423e+00,-5.649988e+00,5.851070e-01,-5.337215e+00,1.672716e+00,-5.293217e+00,5.406546e+00,3.395460e+00,4.991247e+00,1.785692e+00,3.300326e+00,-4.600908e+00,-2.577841e+00,-6.880838e-01,-2.043664e+00,6.061366e+00,-4.934798e+00,-2.098093e+00,7.441896e-02,-5.963899e+00,3.122744e+00,4.313839e+00,-4.315486e+00,-4.787087e-02,5.484700e+00,-4.565164e+00,2.400024e+00,-5.233042e+00,5.626264e+00,-6.945571e+00,1.166851e+00,2.646486e+00,3.424596e+00,5.391860e+00,3.665764e+00};

double Mono_F10_SelfAdaptive::getMaximum (const int i) const {
  if (i == getNumberOfVar() - 1) {
    return 1;
  }
  return 15;
}

double Mono_F10_SelfAdaptive::getMinimum (const int i) const {
  if (i == getNumberOfVar() - 1) {
    return 0;
  }
  return -15;
}

bool Mono_F10_SelfAdaptive::init (const vector<string> &params){

  if (params.size() != 1) {
    cerr << "Error in init: Number of variables must be indicated" << endl;
    exit(-1);
  }
	
  // The total number of variables is equal to the number of variables of the problem plus the
  // the variable which represents the parameter
  setNumberOfVar(atoi(params[0].c_str()) + 1);

	setNumberOfObj(NOBJ);
	return true;
}

// Evaluacion
void Mono_F10_SelfAdaptive::evaluate (void) {
	long double sum = 0.0;
	long double currentGen;
	long double nextGen;

	currentGen = getVar(0)-f10[0];

  // The evaluation only considers the variables of the problem
	for (int i = 1; i < getNumberOfVar() - 1; i++){
		nextGen = getVar(i)-f10[i];
		sum += currentGen * currentGen + 2.0 * nextGen * nextGen;
		sum += -0.3 * cos(3.0 * M_PI * currentGen) -0.4 * cos(4.0 * M_PI * nextGen) + 0.7;
		currentGen = nextGen;
	}
	setObj(0, sum);
}

// Clonacion
Individual* Mono_F10_SelfAdaptive::clone (void) const{
	return new Mono_F10_SelfAdaptive();
}
