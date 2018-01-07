using System;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using SixLabors.ImageSharp;

namespace smallptSharp
{
    struct Vec
    {
        public double x, y, z;          // Also used as RGB color presentation as Rgb(double r, double g, double b)
                                        // RGB implention will be added later.

        public Vec(double _x = 0, double _y = 0, double _z = 0)
        {
            x = _x;
            y = _y;
            z = _z;
        }

        public static Vec operator +(Vec v1, Vec v2) =>
            new Vec(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);

        public static Vec operator -(Vec v1, Vec v2) =>
            new Vec(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);

        public static Vec operator *(Vec v, double b) =>           // Dot product here, not cross product
            new Vec(v.x * b, v.y * b, v.z * b);

        public static Vec operator %(Vec v1, Vec v2) =>
            new Vec(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);

        public Vec Mult(Vec v) =>
            new Vec(x * v.x, y * v.y, z * v.z);

        public Vec Norm() =>
            this * (1 / Math.Sqrt(x * x + y * y + z * z));   // Return the length of this vector.

        public double Dot(Vec v) =>
            x * v.x + y * v.y + z * v.z;
    }

    struct Ray
    {
        public Vec o, d;
        public Ray(Vec _o,Vec _d)
        {
            o = _o;
            d = _d;
        }
    }

    enum Refl_t                                              // Material types, used in radiance()
    {
        DIFF,
        SPEC,
        REFR
    }

    class Sphere
    {
        public double rad;     // radius
        public Vec p, e, c;    // position, emission, color
        public Refl_t refl;
        public Sphere(double _rad, Vec _p, Vec _e, Vec _c, Refl_t _refl)
        {
            rad = _rad;
            p = _p;
            e = _e;
            c = _c;
            refl = _refl;
        }

        public double Intersect(Ray r)      // returns distance, 0 if nohit
        {
            Vec op = p - r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
            double t, eps = 1e-4, b = op.Dot(r.d), det = b * b - op.Dot(op) + rad * rad;
            if (det < 0)
                return 0;
            else
                det = Math.Sqrt(det);
            return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
        }
    }

    class SceneLoader
    {
        public static StreamReader SceneReader;

        public async static Task<Ray> GetCamera()
        {
            string line = await SceneReader.ReadLineAsync();
            var cam = Array.ConvertAll(Regex.Matches(line, @"\.?\-?\d+\.?\d*").Cast<Match>().Select(m => m.Value).ToArray(), double.Parse);
            return new Ray(new Vec(cam[0], cam[1], cam[2]), (new Vec(cam[3], cam[4], cam[5])- new Vec(cam[0], cam[1], cam[2])).Norm());
        }

        public async static Task<int> GetCount()
        {
            return int.Parse((await SceneReader.ReadLineAsync()).Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries)[1]);
        }

        public async static Task<Sphere[]> GetScene()
        {
            string line;
            System.Collections.Generic.List<Sphere> spheres = new System.Collections.Generic.List<Sphere>();
            double[] singleSphere = new double[] { };
            while ((line = await SceneReader.ReadLineAsync()) != null)
            {
                singleSphere = Array.ConvertAll(Regex.Matches(line, @"\ \.?\-?\d+\.?\d*").Cast<Match>().Select(m => m.Value).ToArray(), double.Parse);
                spheres.Add(new Sphere(singleSphere[0], new Vec(singleSphere[1], singleSphere[2], singleSphere[3]), new Vec(singleSphere[4], singleSphere[5], singleSphere[6]), new Vec(singleSphere[7], singleSphere[8], singleSphere[9]), (Refl_t)singleSphere[10]));
            }

            return spheres.ToArray();
        }
    }

    class Program
    {
        static Random rand = new Random();

        public static Sphere[] spheres = new Sphere[]       // Define render scene here.
        {
            new Sphere(1e5,  new Vec(1e5+1, 40.8, 81.6),     new Vec(), new Vec(.75, .25, .25), Refl_t.DIFF),    // Left
            new Sphere(1e5,  new Vec(-1e5+99, 40.8, 81.6),   new Vec(), new Vec(.25, .25, .75), Refl_t.DIFF),    // Right
            new Sphere(1e5,  new Vec(50, 40.8, 1e5),         new Vec(), new Vec(.75, .75, .75), Refl_t.DIFF),    // Back
            new Sphere(1e5,  new Vec(50, 40.8, -1e5+170),    new Vec(), new Vec(),              Refl_t.DIFF),    // Front
            new Sphere(1e5,  new Vec(50, 1e5, 81.6),         new Vec(), new Vec(.75, .75, .75), Refl_t.DIFF),    // Bottom
            new Sphere(1e5,  new Vec(50, -1e5+81.6, 81.6),   new Vec(), new Vec(.75, .75, .75), Refl_t.DIFF),    // Top
            new Sphere(16.5, new Vec(27, 16.5, 47),          new Vec(), new Vec(1, 1, 1)*.999,  Refl_t.SPEC),    // Mirror
            new Sphere(16.5, new Vec(73, 16.5, 78),          new Vec(), new Vec(1, 1, 1)*.999,  Refl_t.REFR),    // Glass
            new Sphere(1.5,  new Vec(50, 81.6 - 16.5, 81.6), new Vec(4, 4, 4)*100, new Vec(),   Refl_t.DIFF)     // Light
                     //Radius    Position                        Emission   Color               ReflectionType
        };

        public static int numSpheres = spheres.Length;

        public static double Clamp(double x) =>
            x < 0 ? 0 : x > 1 ? 1 : x;

        public static int ToInt(double x) =>
            (int)(Math.Pow(Clamp(x), 1 / 2.2) * 255 + .5);

        static bool Intersect(Ray r, ref double t, ref int id)
        {
            double inf = t = 1e20, d;
            for (int i = spheres.Length - 1; i >= 0; i--)
            {
                d = spheres[i].Intersect(r);
                if (d > 0 && d < t)
                {
                    t = d; id = i;
                }
            }
            return t < inf;
        }

        static Vec Radiance(Ray r, int depth)
        {
            double t = 0;               // Distance to intersection
            int id = 0;                 // ID of intersected object
            if (!Intersect(r, ref t, ref id))   // If miss, return black
                return new Vec();
            Sphere obj = spheres[id];   // The hit object
            Vec x = r.o + r.d * t, n = (x - obj.p).Norm(), nl = n.Dot(r.d) < 0 ? n : n * -1, f = obj.c;
            double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z;

            if (depth > 100)
                return obj.e;           // In order to prevent stack overflow.

            if (++depth > 5)
                if (rand.NextDouble() < p)
                    f = f * (1 / p);
                else
                    return obj.e;

            if (obj.refl == Refl_t.DIFF)        // Ideal diffuse reflection
            {
                double r1 = 2 * Math.PI * rand.NextDouble(), r2 = rand.NextDouble(), r2s = Math.Sqrt(r2);
                Vec w = nl, u = ((Math.Abs(w.x) > .1 ? new Vec(0, 1) : new Vec(1)) % w).Norm(), v = w % u;
                Vec d = (u * Math.Cos(r1) * r2s + v * Math.Sin(r1) * r2s + w * Math.Sqrt(1 - r2)).Norm();

                // Loop over any lights
                Vec e = new Vec();
                for (int i = 0; i < numSpheres; i++)
                {
                    Sphere s = spheres[i];
                    if (s.e.x <= 0 && s.e.y <= 0 && s.e.z <= 0)         // Skip non-lights
                        continue;

                    Vec sw = s.p - x, su = ((Math.Abs(sw.x) > .1 ? new Vec(0, 1) : new Vec(1)) % sw).Norm(), sv = sw % su;
                    double cos_a_max = Math.Sqrt(1 - s.rad * s.rad / (x - s.p).Dot(x - s.p));
                    double eps1 = rand.NextDouble(), eps2 = rand.NextDouble();
                    double cos_a = 1 - eps1 + eps1 * cos_a_max;
                    double sin_a = Math.Sqrt(1 - cos_a * cos_a);
                    double phi = 2 * Math.PI * eps2;
                    Vec l = su * Math.Cos(phi) * sin_a + sv * Math.Sin(phi) * sin_a + sw * cos_a;
                    l = l.Norm();
                    if (Intersect(new Ray(x, l), ref t, ref id) && id == i)
                    {
                        double omega = 2 * Math.PI * (1 - cos_a_max);
                        e = e + f.Mult(s.e * l.Dot(nl) * omega) * (1 / Math.PI);    // 1/PI for brdf
                    }
                }

                return obj.e + e + f.Mult(Radiance(new Ray(x, d), depth));
            }
            else if (obj.refl == Refl_t.SPEC)      //Ideal specular reflection
            {
                return obj.e + f.Mult(Radiance(new Ray(x, r.d - n * 2 * n.Dot(r.d)), depth));
            }

            Ray reflRay = new Ray(x, r.d - n * 2 * n.Dot(r.d));        // Ideal dielectric refraction
            bool into = n.Dot(nl) > 0;
            double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.Dot(nl), cos2t;
            if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)         // Total internal reflection
                return obj.e + f.Mult(Radiance(reflRay, depth));
            Vec tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + Math.Sqrt(cos2t)))).Norm();
            double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : tdir.Dot(n));
            double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
            return obj.e + f.Mult(depth > 2 ? (rand.NextDouble() < P ?   // Russian roulette
              Radiance(reflRay, depth) * RP : Radiance(new Ray(x, tdir), depth) * TP) :
              Radiance(reflRay, depth) * Re + Radiance(new Ray(x, tdir), depth) * Tr);
        }

        static async Task Main(string[] args)
        {
            int w = 800, h = 600, samps = 1, sp = 2;
            string scene = "scene.scn";
            string[] arg;
            if (args.Length > 0)
            {
                foreach(string s in args)
                {
                    arg = s.Split(":");
                    switch (arg[0])
                    {
                        case "resolution":
                            w = int.Parse(arg[1].Split("*")[0]);
                            h = int.Parse(arg[1].Split("*")[1]);
                            break;
                        case "samples":
                            samps = int.Parse(arg[1])/4;
                            break;
                        case "scene":
                            scene = arg[1];
                            break;
                        case "subpixels":
                            sp = int.Parse(arg[1]);
                            break;
                    }
                }
            }
            await Console.Out.WriteLineAsync("smallpt C# implementation by net2cn");
            await Console.Out.WriteLineAsync("Powered by .NET Core 2.0");
            Random rand = new Random();
            Ray cam = new Ray(new Vec(50, 52, 295.6), new Vec(0, -0.042612, -1).Norm());
            //using (StreamWriter sw = new StreamWriter("s.scn"))
            //{
            //    await sw.WriteLineAsync($"camera {cam.o.x} {cam.o.y} {cam.o.z} {cam.d.x + cam.o.x} {cam.d.y + cam.o.y} {cam.d.z + cam.o.z}");
            //    await sw.WriteLineAsync($"size {numSpheres}");
            //    foreach (Sphere s in spheres)
            //    {
            //        await sw.WriteLineAsync($"sphere {s.rad} {s.p.x} {s.p.y} {s.p.z} {s.e.x} {s.e.y} {s.e.z} {s.c.x} {s.c.y} {s.c.z} {(int)s.refl}");
            //    }
            //}
            if (scene != null)
            {
                SceneLoader.SceneReader = new StreamReader(scene);
                cam = await SceneLoader.GetCamera();
                numSpheres = await SceneLoader.GetCount();
                spheres = await SceneLoader.GetScene();
            }

            Vec cx = new Vec(w * .5135 / h), cy = (cx % cam.d).Norm() * .5135;      // .5135: Fov
            var c = Enumerable.Repeat(new Vec(), w * h).ToArray();
            for (int y = 0; y < h; y++)     // Loop over image rows
            {
                await Console.Out.WriteAsync($"\rRendering {scene}({w}*{h} image, {samps * 4} spp, {sp} subpixels anti-aliasing) {100 * y / (h - 1)}%");
                for (int x = 0; x < w; x++) // Loop cols
                {
                    for (int sy = 0, i = (h - y - 1) * w + x; sy < sp; sy++)     // subpixel rows
                    {
                        for (int sx = 0; sx < sp; sx++)                          // subpixel cols
                        {
                            Vec r = new Vec();
                            for (int s = 0; s < samps; s++)
                            {
                                double r1 = 2 * rand.NextDouble(), dx = r1 < 1 ? Math.Sqrt(r1) - 1 : 1 - Math.Sqrt(2 - r1);
                                double r2 = 2 * rand.NextDouble(), dy = r2 < 1 ? Math.Sqrt(r2) - 1 : 1 - Math.Sqrt(2 - r2);
                                Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                                        cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
                                d = d.Norm();
                                r = r + Radiance(new Ray(cam.o + d * 140, d), 0) * (1.0 / samps);
                            } // Camera rays are pushed ^^^^^ forward to start in interior
                            c[i] = c[i] + new Vec(Clamp(r.x), Clamp(r.y), Clamp(r.z)) * .25;
                        }
                    }
                }
            }

            await Console.Out.WriteAsync("...Done!\n");
            await Console.Out.WriteLineAsync("Writing image to file...");
            // Write image to PPM file. Will implement write to PNG later.
            //using (StreamWriter sw = new StreamWriter("image.ppm"))
            //{
            //    sw.Write("P3\r\n{0} {1}\r\n{2}\r\n", w, h, 255);
            //    for (int i = 0; i < w * h; i++)
            //        sw.Write("{0} {1} {2}\r\n", ToInt(c[i].x), ToInt(c[i].y), ToInt(c[i].z));
            //    sw.Close();
            //}
            using (Image<Rgba32> image = new Image<Rgba32>(w, h))
            {
                int i = 0;
                for (int y = 0; y < h; y++)
                    for (int x = 0; x<w; x++)
                    {
                        image[x, y] = Rgba32.FromHex(ToInt(c[i].x).ToString("x2") + ToInt(c[i].y).ToString("x2") + ToInt(c[i].z).ToString("x2"));
                        i++;
                    }
                image.Save("image.png");
            }
        }
    }
}
