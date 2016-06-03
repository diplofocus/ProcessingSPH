float H = 5e-2;
float H2 = H*H;
float H8 = (H2 * H2) * (H2 * H2);
float g = 9.81;
float dt = 5e-4;
float scale = 800;

float rho0 = 1000;
float k = 1e3;
float mu = 0.1;

//int nBorder = 400;
int nParticles = 30000;
int n = nParticles; //+ nBorder;


Particle[] particles;

void setup()
{
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
  //      a++;
  //  }
  //}


  NormalizeMass();
  LeapfrogStart();
}

void draw()
{
  //println(random(0, 1));
  background(0);
  ComputeDensity();
  ComputeAcceleration();
  Constrain();
  Leapfrog();

  for (Particle p : particles)
  {
    p.Draw();
  }

  //println(particles[2].ay);
  Flush();
  if (frameCount % 50 == 0)
  {
    saveFrame("#####.png");
  }
}


void Flush()
{
  for (int i = 0; i < n; i++)
  {
    particles[i].rho = 0;
    particles[i].ax = 0;
    particles[i].ay = 0;
  }
}