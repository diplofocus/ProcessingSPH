class Particle //<>//
{
  float x, y, vx, vy, ax, ay, vhx, vhy;
  double[] r;
  float m, h;
  float rho;
  float p;
  int isBorder;

  Particle()
  {
    r = new double[2];
    m = 1;
    h = H;
    x = random(0, 1);
    y = random(0, 1);
    vx = 0;
    vy = 0;
    ax = 0;
    ay = 0;
    vhx = 0;
    vhy = 0;
    isBorder = 0;
    r[0] = x;
    r[1] = y;
  }

  void Draw()
  {
    noStroke();
    fill(30, 30, 255);
    rect(x*scale, y*scale, H*scale/4, H*scale/4);
  }
}

void ComputeDensity()
{
  float C = 4*particles[0].m / 3.1415926 / H8;
  for (int i = 0; i < n; i++)
  {
    particles[i].rho += 4*particles[i].m / 3.1415926 / H2;

    for (Object o : tree.nearest(particles[i].r, 32)) 
    {
      Particle p2 = (Particle) o;

      float dx = particles[i].x - p2.x;
      float dy = particles[i].y - p2.y;
      float dr2 = dx*dx + dy*dy;
      float z = H2 - dr2;
      if (z > 0)
      {
        float rho_ij = C*z*z*z;
        particles[i].rho += rho_ij;

      }
    }

    //for (int j = i+1; j < n; j++)
    //{
    //  float dx = particles[i].x - particles[j].x;
    //  float dy = particles[i].y - particles[j].y;
    //  float dr2 = dx*dx + dy*dy;
    //  float z = H2 - dr2;
    //  if (z > 0)
    //  {
    //    float rho_ij = C*z*z*z;
    //    particles[i].rho += rho_ij;
    //    particles[j].rho += rho_ij;
    //  }
    //}
  }
}

void ComputeAcceleration()
{
  for (int i = 0; i < n; i++)
  {
    particles[i].ax = 0;
    particles[i].ay = g;
  }

  float C0 = particles[0].m / 3.1415926 / (H2*H2);
  float Cp = 15*k;
  float Cv = -40*mu;

  for (int i = 0; i < n; i++)
  {

    for (Object o : tree.nearest(particles[i].r, 32)) 
    {
      Particle p2 = (Particle) o;
      float dx = particles[i].x - p2.x;
      float dy = particles[i].y - p2.y;
      float dr2 = dx*dx + dy*dy;
      if (dr2 < H2)
      {
        float q = sqrt(dr2)/H;
        float u = 1-q;
        float w0 = C0 * u/particles[i].rho/p2.rho;
        float wp = w0 * Cp * (particles[i].rho + p2.rho - 2*rho0) * u/q;
        float wv = w0 * Cv;
        float dvx = particles[i].vx - p2.vx;
        float dvy = particles[i].vy - p2.vy;
        particles[i].ax += (wp*dx + wv*dvx);
        particles[i].ay += (wp*dy + wv*dvy);
      }
    }

    //for (int j = i+1; j < n; j++)
    //{
    //  float dx = particles[i].x - particles[j].x;
    //  float dy = particles[i].y - particles[j].y;
    //  float dr2 = dx*dx + dy*dy;
    //  if (dr2 < H2)
    //  {
    //    float q = sqrt(dr2)/H;
    //    float u = 1-q;
    //    float w0 = C0 * u/particles[i].rho/particles[j].rho;
    //    float wp = w0 * Cp * (particles[i].rho + particles[j].rho - 2*rho0) * u/q;
    //    float wv = w0 * Cv;
    //    float dvx = particles[i].vx - particles[j].vx;
    //    float dvy = particles[i].vy - particles[j].vy;
    //    particles[i].ax += (wp*dx + wv*dvx);
    //    particles[i].ay += (wp*dy + wv*dvy);
    //    particles[j].ax -= (wp*dx + wv*dvx);
    //    particles[j].ay -= (wp*dy + wv*dvy);
    //  }
    //}
  }
}

void Integrate()
{
  for (int i = 0; i < n; i++)
  {
    if (particles[i].isBorder == 0)
    {
      particles[i].vx += particles[i].ax*dt;
      particles[i].vy += particles[i].ay*dt;

      tree.delete(particles[i].r);

      particles[i].x += particles[i].vx*dt;
      particles[i].y += particles[i].vy*dt;

      particles[i].r[0] = particles[i].x;
      particles[i].r[1] = particles[i].y;
      tree.insert(particles[i].r, particles[i]);
    }
  }
}

void Constrain()
{
  for (int i = 0; i < n; i++)
  {
    if (particles[i].x < XMIN)
    {
      particles[i].x = 0;
      particles[i].vx *= -D;
      particles[i].vhx *= -D;
      particles[i].x -= particles[i].vx * (1-D) * dt;
    }

    if (particles[i].x > XMAX)
    {
      particles[i].x = XMAX;
      particles[i].vx *= -D;
      particles[i].vhx *= -D;
      particles[i].x -= particles[i].vx * (1-D) * dt;
    }

    if (particles[i].y < YMIN)
    {
      particles[i].y = 0;
      particles[i].vy *= -D;
      particles[i].vhy *= -D;
      particles[i].y -= particles[i].vy * (1-D) * dt;
    }

    if (particles[i].y > YMAX)
    {
      particles[i].y = YMAX;
      //println("aaaaa");
      particles[i].vy *= -D;
      particles[i].vhy *= -D;
      particles[i].y -= particles[i].vy * (1-D) * dt;
    }
  }
}

void NormalizeMass()
{
  float m = 1;
  ComputeDensity();
  float rho2s = 0;
  float rhos = 0;
  for (int i = 0; i < n; i++)
  {
    rho2s += particles[i].rho * particles[i].rho;
    rhos += particles[i].rho;
  }
  m *= rho0*rhos/rho2s;
  for (int i = 0; i < n; i++)
  {
    particles[i].m = m;
  }
}


void LeapfrogStart()
{
  for (int i = 0; i < n; i++)
  {
    particles[i].vhx = particles[i].vx + particles[i].ax * dt /2.0;
    particles[i].vhy = particles[i].vy + particles[i].ay * dt /2.0;

    particles[i].vx += particles[i].ax * dt;
    particles[i].vy += particles[i].ay * dt;

    particles[i].x += particles[i].vhx * dt;
    particles[i].y += particles[i].vhy * dt;

    Constrain();
  }
}

void Leapfrog()
{
  for (int i = 0; i < n; i++)
  {
    particles[i].vhx += particles[i].ax * dt;
    particles[i].vhy += particles[i].ay * dt;

    particles[i].vx = particles[i].vhx + particles[i].ax * dt/2.0;
    particles[i].vy = particles[i].vhy + particles[i].ay * dt/2.0;

    particles[i].x += particles[i].vx * dt;
    particles[i].y += particles[i].vy * dt;

    Constrain();
  }
}

//void Reflect(int w, float barrier, int c)
//{
//  if (particles[w].vx == 0 && particles[w].vy == 0)
//    return;

//  switch(c)
//  {
//  case 0:
//    {
//      particles[w].vx *= -0.5;
//      particles[w].vhx *= -0.5;
//    }
//    break;

//  case 1:
//    {
//      particles[w].vy *= -0.5;
//      particles[w].vhy *= -0.5;
//    }
//    break;

//  case 2:
//    {
//      particles[w].vx *= -0.5;
//      particles[w].vhx *= -0.5;
//    }
//    break;

//  case 3:
//    {
//      particles[w].vy *= -0.5;
//      particles[w].vhy *= -0.5;
//    }
//    break;
//  }
//}

//void Constrain()
//{
//  float damp = 0.25;
//  for (int i = 0; i < n; i++)
//  {
//    if (particles[i].x > 2)
//    {
//      Reflect(i, 0, 2);
//    }

//    if (particles[i].x < 0)
//    {
//      Reflect(i, 0, 0);
//    }

//    if (particles[i].y > 1)
//    {
//      Reflect(i, 4, 1);
//    }

//    if (particles[i].y < 0)
//    {
//      Reflect(i, 2, 3);
//    }
//  }
//}