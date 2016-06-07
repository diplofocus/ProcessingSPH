import net.sf.javaml.core.kdtree.KDTree;

KDTree tree;

float H = 5e-2;
float H2 = H*H;
float H8 = (H2 * H2) * (H2 * H2);
float g = 9.81;
float dt = 5e-4;
float scale = 800;
float D = 0.95;
int m = 16;

float rho0 = 1000;
float k = 1e3;
float mu = 10;

//int nBorder = 400;
int nParticles = 3000;
int n = nParticles; //+ nBorder;

int XMIN = 0;
int XMAX = 2;
int YMIN = 0;
int YMAX = 1;


Particle[] particles;

void setup()
{
  tree = new KDTree(2);
  rectMode(CENTER);
  size(1600, 800);
  particles = new Particle[n];
  for (int i = 0; i < n; i++)
  {
    particles[i] = new Particle();
  }

  //for (int i = nParticles; i < n; i++)
  //{
  //  particles[i].isBorder = 1;
  //}


  //int a = 0;
  //for (int i = 0; i < 40; i++)
  //{
  //  for (int j = 0; j < 40; j++)
  //  {

  //    particles[a].x = i*H;
  //    particles[a].y = j*H;
  //    a++;
  //  }
  //}

  for (Particle p : particles)
  {
    p.r[0] = p.x;
    p.r[1] = p.y;
    //tree.insert(p.r, p);
  }
  MakeTree();
  NormalizeMass();
  PurgeTree();
  LeapfrogStart();
  MakeTree();
}

void draw()
{

  //println(random(0, 1));
  background(0);
  ComputeDensity();
  ComputeAcceleration();
  //PurgeTree();
  Constrain();
  Leapfrog();
  for (Particle p : particles)
  {
    p.r[0] = p.x;
    p.r[1] = p.y;
  }
  //MakeTree();

  for (Particle p : particles)
  {
    p.Draw();
  }
  //println(particles[5].r[0]);

  //println(particles[2].ay);
  Flush();
  if (frameCount % 50 == 0)
  {
    saveFrame("#####.png");
  }
}


void Flush()
{
  tree = new KDTree(2);
  for (int i = 0; i < n; i++)
  {
    particles[i].rho = 0;
    particles[i].ax = 0;
    particles[i].ay = 0;
    tree.insert(particles[i].r, particles[i]);
    particles[i].c = 0;
  }
}
