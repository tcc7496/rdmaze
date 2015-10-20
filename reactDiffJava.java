import java.util.*;
import java.util.Random;
import java.applet.*;
import java.awt.*;
import java.awt.event.*;
import java.math.*;
import java.lang.*;

/*      <applet code="reactDiffJava" width=550 height=550>
         </applet>
*/

public class maze_generator2 extends Applet implements Runnable, ActionListener { 
	int nx, ny;
	double eps, dt, dx, diff, delta, a, b, gamma, dlap, nf, t, beta, ntime;
	double[][] u= new double[101][101];
	double[][] v= new double[101][101];
	double[][] ut= new double[101][101];
	double[][] phi= new double[101][101];
	int xscale=200, yscale=200;
public void init(){
    nf = 19;
    eps = 0.08;
    t = 0;
    gamma = 0.8;
    beta = 0.7;
    dt = 0.005;
    dx = 0.25;
    diff = 1.0;
    dlap = diff*dt/(dx*dx);
    ntime = 15000;
    initGameElements();
 }
public void initGameElements(){
	u = new double[101][101];
	v = new double[101][101];
	ut = new double[101][101];
	phi = new double[101][101];
	nf = 19;
    eps = 0.08;
    t = 0;
    gamma = 0.8;
    beta = 0.7;
    dt = 0.005;
    dx = 0.25;
    diff = 1.0;
    dlap = diff*dt/(dx*dx);
    ntime = 15000;

		for(l = 0; l < 101; l++){
			for(m = 0; m < 101; m++){
				u[l][m] = -1.199;
				ut[l][m] = u[l][m];
				v[l][m] = -0.6242;
				phi[l][m] = 1;
			}
		}

		for(p = 47; p < 63; p++){
			for(q = 47; q < 63; q++){
				u[p][q] = 1;
			}
		}
}





}
