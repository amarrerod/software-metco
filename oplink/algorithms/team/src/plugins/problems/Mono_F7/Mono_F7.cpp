/***********************************************************************************
 * AUTORES
 *   - Eduardo Manuel Segredo Gonz�lez
 *   - Carlos Segura Gonz�lez
 * 
 * FECHA
 *    Febrero 2011
 ************************************************************************************/

#include <math.h>
#include "Mono_F7.h"

double f7[1000] = {3.699709e+00,4.351871e+00,-3.963818e+00,-2.493211e+00,5.037569e-01,2.681006e+00,-2.579118e+00,3.736011e+00,4.680081e+00,-2.016519e+00,4.648024e+00,4.410581e+00,2.538111e+00,1.618088e+00,1.104366e+00,3.289030e-01,4.946150e+00,-4.046462e+00,-2.364280e+00,-4.754897e+00,8.307132e-02,-1.190161e+00,-2.496306e+00,3.512223e+00,1.543734e+00,3.084913e+00,-4.661524e+00,5.886419e-01,-6.331781e-01,4.834074e-01,4.809783e+00,-2.014026e+00,1.966675e+00,-1.966087e+00,-3.686643e+00,1.702272e+00,-4.322699e+00,8.937012e-01,-1.728210e+00,8.429327e-01,1.274640e+00,-3.605217e+00,2.417219e+00,-3.403977e+00,-1.403683e+00,-1.436441e+00,-3.278730e+00,-3.182312e-01,3.797633e+00,1.769595e+00,
2.543461e+00,3.364133e+00,-1.775673e+00,-4.887483e+00,-4.370645e-01,7.881307e-01,8.115332e-01,-3.712902e+00,7.506810e-01,2.191919e+00,2.639957e+00,-2.723791e+00,2.515829e+00,-2.719982e+00,3.274094e+00,4.239466e+00,-8.353419e-01,5.937957e-01,2.941319e+00,4.801043e+00,4.972368e+00,-4.601751e+00,-2.104079e+00,-6.322118e-01,4.954399e+00,4.043619e+00,-6.701446e-03,4.891432e+00,1.754707e+00,1.734330e+00,4.789139e+00,2.548223e+00,-3.791961e+00,1.386961e+00,2.386926e+00,1.785687e+00,-9.228598e-01,-1.029663e+00,-3.673604e+00,1.613214e+00,2.518770e+00,-4.743181e-01,-2.552215e+00,-3.585540e+00,-9.702251e-01,3.892420e+00,4.379854e+00,1.041770e-01,-1.913105e+00,4.446294e+00,
2.478954e+00,-9.291761e-01,-4.977900e+00,4.483526e-01,3.554588e+00,-1.232989e+00,3.811400e+00,-8.604531e-01,3.735492e+00,-7.340429e-01,-1.761808e+00,1.265607e+00,-4.332418e+00,1.481440e+00,-2.345555e+00,-4.472322e-01,6.436119e-01,-2.040895e-01,9.162494e-01,7.598824e-01,7.632996e-01,1.789720e+00,-4.731494e+00,-3.806022e+00,-4.240923e+00,-2.090942e+00,-1.244739e+00,2.375820e+00,-1.831512e+00,1.159392e+00,-2.246483e+00,1.924547e+00,3.511957e+00,-1.533483e+00,-9.579426e-01,-1.869493e+00,-2.795329e+00,-4.044361e+00,1.121844e+00,-1.306600e+00,4.934400e+00,4.947103e+00,4.581485e+00,4.173404e+00,-2.040958e+00,4.025699e-01,1.406554e+00,-4.803102e+00,-4.067694e+00,-3.293897e+00,
-1.655845e+00,-1.440755e+00,-1.631603e+00,9.116558e-01,3.328525e-01,-1.582515e+00,1.603256e+00,-2.216722e+00,2.297742e+00,3.293911e+00,-2.985253e+00,9.572424e-01,-4.857932e+00,-1.142152e+00,-2.840986e+00,4.706327e+00,-3.513329e+00,-9.652393e-01,-2.838829e+00,4.758874e+00,-2.264914e+00,-5.652813e-01,-2.141403e+00,-3.954953e+00,-5.809671e-01,-3.200707e+00,2.409460e+00,1.559364e+00,2.850187e+00,2.355555e+00,1.938146e+00,-4.171654e+00,8.861805e-01,2.133378e+00,-2.195231e-01,-3.387101e+00,-3.638647e+00,7.158923e-01,1.871076e+00,2.306677e+00,4.938924e+00,-4.577942e+00,1.197990e+00,9.103813e-01,4.217289e-01,3.017394e+00,-1.673030e+00,3.292847e+00,3.683001e+00,-3.616799e+00,
-8.213508e-01,-3.226428e-01,-4.771577e+00,-8.501453e-01,2.963055e+00,3.798880e+00,2.483576e+00,9.085046e-01,3.676223e+00,3.341080e+00,-4.065719e+00,-3.066258e+00,-1.563804e+00,-3.865908e+00,-4.185267e+00,3.901124e-02,-2.802633e-01,4.645685e+00,3.591044e+00,-4.145506e+00,3.728630e+00,4.714492e+00,4.583040e+00,3.644515e+00,1.021904e+00,-3.621365e+00,1.029173e+00,8.054410e-01,-4.711201e+00,3.118483e+00,4.690942e-01,4.822737e+00,4.273555e+00,-7.418578e-01,-3.434397e+00,-4.083407e+00,3.034074e+00,-3.707426e+00,-2.276230e+00,-2.454403e+00,3.485837e+00,-1.365673e+00,-3.610931e+00,-2.654396e+00,-8.053430e-01,-3.038668e+00,-1.907013e+00,-7.274325e-01,-6.207556e-01,4.871475e+00,
-3.186744e+00,1.548188e+00,4.621232e+00,8.886594e-01,1.477014e+00,-1.126102e+00,2.655670e+00,-2.755373e+00,3.440405e+00,-2.705599e+00,2.696215e+00,-9.231469e-02,-1.290060e+00,1.427450e+00,-4.059305e+00,-1.990343e+00,-6.797871e-01,2.692098e+00,3.761556e+00,1.658898e+00,-3.717188e+00,8.453135e-01,-2.102973e+00,-1.137698e+00,2.103071e+00,1.438892e+00,3.540443e+00,-4.169609e+00,1.176913e+00,-2.350457e+00,-3.924295e+00,4.866587e+00,-5.334477e-01,4.398915e+00,-1.475992e+00,2.434978e+00,-4.140821e-01,-3.474255e+00,3.241112e+00,-5.291902e-01,-3.053990e+00,-1.695970e+00,-3.617261e+00,-2.990056e+00,2.545772e+00,-2.094696e+00,-2.228949e+00,1.811512e+00,-6.722944e-01,4.312503e+00,
1.519814e+00,3.815882e+00,2.769014e+00,4.037513e+00,4.855495e+00,4.904709e+00,2.022051e+00,3.817240e+00,3.614208e+00,-2.376506e+00,4.861251e+00,-3.132166e+00,-4.940240e+00,3.385281e+00,-1.390687e+00,1.695774e+00,4.657617e+00,-3.922482e-01,-1.484143e+00,3.351780e+00,4.303722e+00,2.200700e+00,3.194937e+00,-1.230398e+00,-2.976247e+00,4.711201e+00,3.168566e+00,-1.570358e+00,4.181576e-01,-1.166409e+00,2.369476e+00,-1.121382e+00,4.149792e+00,-4.478954e+00,3.567921e+00,4.138181e+00,1.327860e+00,4.012486e+00,-2.782024e+00,1.011344e+00,1.208480e+00,2.054585e+00,-6.217079e-01,4.518560e+00,-4.626274e+00,6.966633e-01,-4.851966e+00,1.307132e+00,3.298253e+00,-2.890998e+00,
-3.837870e+00,-1.683463e+00,-1.254879e+00,2.003942e+00,-4.695011e+00,4.452666e+00,-4.636721e+00,6.819159e-01,-1.389664e+00,-4.002724e+00,3.358713e+00,1.583845e+00,3.269766e+00,-1.833808e+00,4.638052e+00,-3.085347e+00,4.009012e+00,2.526627e+00,2.613221e+00,-4.140576e+00,4.400056e-01,2.201344e+00,-9.580407e-01,4.472217e+00,1.235454e+00,3.378796e+00,-1.015938e+00,4.158573e+00,4.528686e+00,-3.533259e+00,-2.257274e-01,1.397367e+00,1.463457e+00,-1.525038e+00,-4.081195e+00,5.583208e-01,3.023206e+00,4.775792e+00,2.185771e+00,-4.209327e+00,-2.497875e+00,1.638885e+00,-6.077588e-01,-1.076608e+00,1.841511e+00,2.030706e+00,2.548909e+00,1.085977e+00,-1.839130e+00,1.198088e+00,
-1.943412e+00,-4.884920e+00,-1.168117e+00,2.440860e+00,-4.354756e+00,6.940583e-01,-3.538567e+00,9.895382e-01,-1.012983e+00,-6.762158e-01,1.057897e+00,-4.131683e+00,-4.311467e+00,2.618333e+00,-9.692448e-01,-2.828017e+00,-1.514002e+00,6.358531e-01,1.280746e+00,-3.227681e+00,-1.576352e+00,4.595056e+00,3.701838e+00,4.734659e+00,4.309478e+00,6.456287e-01,3.959980e+00,-2.698092e+00,-2.564679e+00,1.957249e+00,-1.622751e+00,2.324688e+00,-4.271356e+00,2.339554e-02,-4.983656e+00,3.424005e+00,2.411953e+00,-3.324344e+00,-1.395462e+00,1.738265e+00,3.002955e+00,-3.362228e+00,-8.388992e-01,1.430692e-01,3.178775e+00,2.876139e+00,-2.842127e-01,-5.609257e-01,-3.610581e+00,2.948461e+00,
-4.962648e+00,2.710346e+00,-2.007598e+00,-1.899170e+00,-4.334071e+00,4.445328e+00,4.041686e+00,9.141627e-01,1.302216e+00,1.464900e+00,3.911803e+00,-1.701082e+00,-4.272476e+00,-3.941214e+00,-2.417429e+00,-1.309408e-01,-4.601330e+00,4.743300e+00,6.128847e-01,-3.432023e-01,-1.503190e+00,-4.127313e+00,1.220595e+00,-3.553566e+00,-2.636343e+00,3.435349e+00,-5.537691e-01,3.756136e+00,-1.564056e+00,2.594090e+00,-5.693568e-01,-1.396625e+00,2.406729e+00,6.645916e-01,-3.896341e+00,-1.176128e+00,-4.256259e+00,3.383390e+00,-3.949519e+00,3.910584e+00,2.340387e+00,3.745142e+00,-4.112706e+00,6.399426e-01,2.977158e+00,-1.745240e+00,4.082189e+00,2.617450e+00,-4.925871e+00,3.115360e+00,
-5.447218e-01,2.070845e+00,4.803746e+00,-1.751906e+00,3.175589e-01,1.177907e+00,3.124856e+00,3.323910e+00,-1.497742e+00,-8.433668e-01,2.724344e-01,-1.837401e-01,-1.177627e+00,-1.028542e+00,3.234151e+00,4.443339e+00,4.418816e+00,-2.957466e+00,4.498589e+00,2.550422e+00,4.719688e+00,-2.577606e+00,-2.001821e-01,-3.019803e+00,1.453766e+00,1.116565e+00,-3.875712e+00,-3.743966e+00,-5.819194e-01,-1.626491e+00,2.840090e+00,-3.285956e+00,4.945310e+00,-1.425867e+00,-4.673037e+00,8.567417e-01,2.211316e+00,3.860110e+00,-2.139806e+00,-3.969434e+00,4.232702e+00,4.020693e+00,3.215693e+00,-7.174749e-01,4.964119e+00,-2.437415e+00,-9.284759e-01,-6.281923e-01,4.433528e-01,-2.205875e-01,
-6.428696e-01,1.991842e+00,3.664515e+00,-4.282784e+00,4.009741e+00,1.382634e+00,4.907762e+00,-7.620111e-01,-2.655838e+00,-2.292140e+00,-1.153062e+00,2.062008e+00,4.000049e+00,-1.396240e-01,1.345121e-01,-7.413396e-01,-1.313855e+00,3.431273e+00,-1.128007e+00,-3.743433e+00,2.103463e+00,3.583257e+00,1.484199e+00,-2.707741e+00,2.387893e+00,-1.208928e+00,2.636357e+00,2.695123e+00,3.925122e+00,2.856280e+00,6.297329e-01,-3.136914e+00,3.824159e+00,-4.857218e+00,2.846728e+00,-1.228171e+00,-2.118000e+00,2.153587e+00,2.377767e+00,-3.264445e+00,-2.671706e+00,-1.940247e+00,-1.349147e+00,-4.202955e+00,2.124274e+00,-2.573684e+00,-3.922223e+00,-2.404839e+00,5.219775e-01,-2.350275e+00,
-3.835685e+00,-1.906873e+00,4.177816e+00,2.462862e+00,3.064059e+00,2.891124e+00,-2.632954e+00,2.798284e+00,4.121599e+00,8.986870e-01,-2.505164e-01,3.340716e+00,-3.611043e+00,-2.005021e+00,-4.303281e-01,2.529008e+00,-3.542376e+00,-1.536900e+00,4.428437e+00,-2.785988e+00,-9.974511e-01,3.815644e+00,2.718217e+00,-1.170582e+00,4.521431e+00,2.957060e+00,2.088729e+00,-4.327236e+00,3.607227e-01,3.067253e+00,-1.420041e+00,-2.395203e+00,4.281510e+00,4.313722e+00,-9.654774e-01,4.520409e+00,4.321774e+00,-1.728910e+00,-7.666888e-01,4.164483e+00,3.993229e+00,1.842127e+00,3.402409e+00,-4.123420e+00,-3.170022e+00,2.941991e+00,6.571128e-01,-5.687546e-01,5.317671e-01,-2.432527e+00,
-3.793992e+00,-4.344463e+00,-2.573894e+00,-4.135079e-01,1.232387e+00,1.274780e+00,-4.971738e+00,4.206596e+00,-2.527678e+00,3.957039e+00,-1.902153e+00,2.577130e+00,5.536151e-01,-4.099835e+00,8.415882e-01,-1.059557e-01,3.459214e+00,-3.804300e+00,1.229922e+00,-4.280459e+00,-2.813438e+00,-4.673989e+00,-4.734421e+00,4.349358e-01,-2.030580e+00,-3.227135e+00,-3.277329e+00,1.929491e+00,3.507069e+00,1.238031e+00,4.919176e+00,1.597052e+00,-6.578411e-01,2.444207e+00,4.305613e+00,-2.213109e+00,-2.598599e+00,-9.330696e-01,-4.577914e+00,-7.650502e-01,1.526067e-01,2.399013e+00,1.460334e+00,9.670180e-01,1.863387e+00,-6.684150e-01,9.493715e-01,-2.650937e+00,-8.440391e-01,-3.641322e+00,
-1.299919e+00,-7.149259e-01,-1.869395e+00,-2.995210e+00,-1.046721e+00,2.444235e+00,3.954714e+00,2.540184e+00,-4.343987e+00,-3.291684e+00,-1.708084e+00,4.006617e+00,1.158286e+00,-1.661433e+00,2.733175e+00,4.155842e+00,-4.426029e+00,9.695249e-01,4.775652e+00,-1.772354e+00,-4.067155e-01,2.428928e+00,2.896110e+00,-2.663219e+00,6.540037e-01,1.636406e+00,-1.184853e+00,-3.625762e+00,-4.082693e+00,-4.008312e+00,-1.086860e+00,-2.487819e+00,-3.261994e+00,3.135051e+00,-6.434649e-02,-1.214313e-01,3.298211e+00,2.175393e+00,3.405910e+00,4.058366e+00,1.539911e+00,1.657274e+00,2.884626e+00,-9.670880e-01,-2.216904e+00,-1.025517e+00,2.873926e+00,-2.498827e+00,-1.639529e+00,6.402437e-02,
3.519198e+00,3.338560e+00,8.723084e-02,1.641896e+00,2.869977e+00,-1.575316e+00,-4.871293e+00,1.978432e-01,1.352733e+00,-3.613116e+00,-1.247610e+00,4.980148e-01,3.227583e+00,-4.209201e+00,5.243724e-01,1.957263e+00,3.869857e+00,3.199587e+00,7.350233e-01,3.755674e+00,-2.230419e+00,2.400343e+00,4.337012e+00,-7.468156e-01,-2.072918e+00,2.950324e+00,2.814957e-01,1.173369e+00,6.317356e-01,1.359966e-01,-4.976121e+00,5.721088e-03,3.415518e+00,2.675978e+00,-9.961486e-01,-4.018074e+00,-2.302013e+00,3.095473e+00,3.994279e+00,2.498799e+00,-1.096691e+00,2.563139e+00,-2.251147e+00,-3.737565e+00,3.388404e+00,2.403102e+00,-4.470873e+00,-4.566738e+00,-1.177879e+00,-2.774000e+00,
6.722384e-01,-2.770701e-01,4.711187e+00,1.859060e+00,-4.066167e+00,-3.053331e+00,1.944400e-01,3.978677e+00,1.536662e+00,4.727355e-01,4.874066e+00,1.180176e+00,-2.719198e+00,-3.928973e+00,2.261384e+00,-1.240916e+00,4.703400e+00,2.020237e-02,2.141501e+00,-3.729456e+00,-1.775267e+00,-3.411428e+00,-1.864452e+00,-4.492245e+00,-4.769826e+00,9.824236e-01,-1.900361e+00,4.947271e+00,-9.629985e-01,3.715437e+00,-1.364581e+00,1.543293e-01,8.645426e-01,1.758979e+00,4.705655e+00,-2.682189e-01,3.245173e+00,-4.529568e+00,-4.914933e+00,-1.921648e+00,1.423220e+00,1.471622e+00,3.601478e+00,-3.891761e+00,1.394860e+00,-3.268464e+00,-2.483982e+00,-4.768425e+00,3.314807e+00,5.791114e-03,
7.262001e-01,2.727769e+00,-2.492700e-01,2.500046e+00,-3.806652e+00,-1.736599e+00,4.627982e+00,-7.969819e-01,4.822870e-01,-4.323371e+00,2.046868e+00,4.729365e+00,2.858450e+00,4.009734e-01,-8.961661e-01,2.892966e-01,-9.504569e-02,2.469234e+00,-4.536165e+00,5.144428e-01,2.033563e+00,3.316894e+00,-3.896061e+00,3.776359e+00,-4.061573e+00,5.088127e-01,-3.019635e+00,7.069991e-01,-4.747460e+00,4.318441e+00,-9.912118e-02,4.864655e+00,1.388390e+00,1.986296e+00,-4.912580e+00,-1.579567e-01,2.626904e+00,-2.170113e+00,-3.036875e+00,-2.655880e+00,-2.024684e+00,-1.991030e+00,1.344428e+00,4.890060e+00,-2.080690e+00,4.037037e+00,-4.059732e-01,-2.897665e+00,2.131025e+00,-8.974966e-01,
1.862617e+00,4.709660e+00,-2.010035e+00,1.705101e+00,-2.799601e+00,-1.650355e+00,3.444018e+00,2.098533e+00,2.291579e+00,2.492147e+00,2.231399e+00,3.324835e+00,-2.589286e+00,-1.258660e+00,1.385519e+00,-3.718322e+00,1.132881e+00,3.653633e+00,2.391408e+00,3.709247e+00,2.292364e+00,-1.227065e+00,9.186863e-01,-1.050699e+00,-3.389580e+00,-1.063079e+00,4.303932e+00,1.803935e+00,-3.568061e+00,3.935506e-01,-2.470901e+00,-2.688722e+00,4.697602e+00,-4.427177e+00,3.944449e+00,-3.471188e+00,3.351402e+00,2.019656e+00,-4.036448e+00,-2.585547e+00,4.695319e+00,2.904135e+00,1.405868e+00,-4.379882e+00,4.960618e+00,-7.224187e-01,1.713686e+00,-3.956899e+00,1.484192e-01,3.592080e+00,
2.116207e+00,3.139407e+00,1.599139e+00,-2.684570e-01,1.756374e+00,4.568376e+00,3.682777e+00,3.692665e+00,4.271328e+00,4.928224e+00,-4.227373e-01,-8.831133e-01,-4.456980e+00,4.161640e+00,4.233675e-01,3.837128e+00,-1.461924e-01,1.289023e+00,2.038731e+00,-8.356500e-01,-4.813382e+00,2.406197e+00,-2.147453e+00,4.377403e+00,4.970204e-01,-4.743608e+00,-2.192759e+00,-2.002766e+00,-4.104548e-01,6.559224e-01,2.689423e+00,2.083099e+00,1.483372e+00,-1.152586e+00,1.530920e+00,2.553853e+00,3.929183e+00,1.913455e+00,-3.064998e+00,2.151388e+00,-4.170792e-01,4.194972e+00,3.942460e+00,-4.477644e-01,-2.083197e+00,-1.827800e+00,-2.399587e+00,-3.944421e+00,-2.484542e+00,-1.370547e+00};


bool Mono_F7::init (const vector<string> &params){

  if (params.size() != 1) {
    cerr << "Error in init: Number of variables must be indicated" << endl;
    exit(-1);
  }
	
  setNumberOfVar(atoi(params[0].c_str()));
	setNumberOfObj(NOBJ);
	return true;
}

// Evaluacion
void Mono_F7::evaluate (void) {
	long double sum, currentGen, prod;
	sum = 0.0;
	prod = 1.0;

	for (int i = 0; i < getNumberOfVar(); i++){   
		currentGen = fabs(getVar(i) - f7[i]);
		sum += currentGen;
		prod *= currentGen;
	}   
	setObj(0, sum + prod);
}

// Clonacion
Individual* Mono_F7::clone (void) const{
	return new Mono_F7();
}
