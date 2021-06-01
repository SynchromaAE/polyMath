/* lcmgcd - find the lowest common multiple and greatest common factor of a list of numbers
 * (c) 2015 Ed Kelly
 */

 #include "m_pd.h"

static t_class *lcmgcd_class;

typedef struct _lcmgcd
{
  t_object x_obj;
  int lcm, gcd, gcdAll, result, b;
  int iterations;
  t_outlet *lcmOut, *gcdOut;
} t_lcmgcd;

static void lcmgcd_calculate_lcm(t_lcmgcd *z) {
  int a, b, t;
 
  a = z->result;
  b = z->b;
 
  while (b != 0) {
    t = b;
    b = a % b;
    a = t;
  }
 
  z->gcd = a;
  z->lcm = (z->result*z->b)/a;
}

static void lcmgcd_calculate_gcd(t_lcmgcd *z) {
  int a, b, t;

  a = z->gcdAll;
  b = z->b;

  while (b != 0) {
    t = b;
    b = a % b;
    a = t;
  }

  z->gcdAll = a;
}

static void lcmgcd_list(t_lcmgcd *z, t_symbol *s, int argc, t_atom *argv)
{
  int i;
  z->iterations = 0;
  z->gcdAll = 1;
  if(argc > 1)
  {
    z->result = (int)atom_getfloat(argv++);
    argc--;
    // make sure the first value is > 0
    while(z->result < 1 && argc)
      {
	z->result = (int)atom_getfloat(argv++);
	argc--;
      }
    for(i=0;i<argc;i++)
      {
	z->b = (int)atom_getfloat(argv++);
	if(z->b > 0)
	  {
	    lcmgcd_calculate_lcm(z);
	    z->result = z->lcm;
	    if(!z->iterations) z->gcdAll = z->gcd;
	    else
	      {
		lcmgcd_calculate_gcd(z);
	      }
	    z->iterations++;
	  }
      }
    outlet_float(z->gcdOut,(float)z->gcdAll);
    outlet_float(z->lcmOut,(float)z->lcm);
 
  }
  else
    {
      post("You need to provide a list");
      post("with more than one number!");
    }
}

static void *lcmgcd_new()
{
  t_lcmgcd *z = (t_lcmgcd *)pd_new(lcmgcd_class);
  z->result = 1;  
  z->lcm = 1;
  z->gcd = 1;
  z->gcdAll = 1;
  z->lcmOut = outlet_new(&z->x_obj, &s_float);
  z->gcdOut = outlet_new(&z->x_obj, &s_float);
  return(void *)z;
}

void lcmgcd_setup(void)
{
  lcmgcd_class = class_new(gensym("lcmgcd"), (t_newmethod)lcmgcd_new, 
			      0, sizeof(t_lcmgcd), 0 ,A_GIMME, 0);
  class_addlist    (lcmgcd_class, lcmgcd_list);
}
