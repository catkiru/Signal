using System.Drawing.Imaging;
using SixLabors.ImageSharp;
using SixLabors.ImageSharp.PixelFormats;
using Color = System.Drawing.Color;

namespace MatrixConvolution;

public static class Tools
{
    static List<Rgb24> colorList = [];

    static Tools()
    {
        var rnd = new Random();
        for (int i = 0; i < 200; i++)
        {
            colorList.Add(new Rgb24((byte)rnd.Next(100, 200), (byte)rnd.Next(100, 200), (byte)rnd.Next(100, 200)));
        }
    }

    public static void SaveMatrixAsBmp(int[,] matrix, List<(int[,] signal, int errors)> signals, string path)
    {
        int height = matrix.GetLength(0) + signals.Sum(x => x.signal.GetLength(0));
        int width = matrix.GetLength(1);

        using var image = new Image<Rgb24>(width, height);

        for (int y = 0; y < matrix.GetLength(0); y++)
        {
            for (int x = 0; x < width; x++)
            {
                int v = matrix[y, x];
                image[x, y] = v == 1 ? new Rgb24(0, 0, 0) : new Rgb24(255, 255, 255);
            }
        }

        for (var i = 0; i < signals.Count; i++)
        {
            var signal = signals[i].signal;
            for (int y = 0; y < signal.GetLength(0); y++)
            {
                for (int x = 0; x < width; x++)
                {
                    int v = signal[y, x];
                    image[x, y + signal.GetLength(0) * i + signal.GetLength(0)] =
                        v == 1 ? colorList[i] : new Rgb24(255, 255, 255);
                }
            }
        }


        image.SaveAsBmp(path);
    }

    public static void SaveMatrixAsBmp2(int[,] matrix, List<(int[,] signal, int errors)> signals, string path)
    {
        int height = matrix.GetLength(0);
        int width = matrix.GetLength(1);

        using var image = new Image<Rgb24>(width, height);

        for (int y = 0; y < matrix.GetLength(0); y++)
        {
            for (int x = 0; x < width; x++)
            {
                int v = matrix[y, x];
                image[x, y] = v == 1 ? new Rgb24(0, 0, 0) : new Rgb24(255, 255, 255);
            }
        }

        for (var i = 0; i < signals.Count; i++)
        {
            var signal = signals[i].signal;
            for (int y = 0; y < signal.GetLength(0); y++)
            {
                for (int x = 0; x < width; x++)
                {
                    int v = signal[y, x];
                    image[x, y] =
                        v == 1 ? colorList[i] : image[x, y];
                }
            }
        }


        image.SaveAsBmp(path);
    }

    public static void SaveMatrixAsBmp(int[,] matrix, string path)
    {
        int height = matrix.GetLength(0);
        int width = matrix.GetLength(1);

        using var image = new Image<L8>(width, height);

        for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
        {
            int v = matrix[y, x];
            image[x, y] = new L8(v == 1 ? (byte)0 : (byte)255);
        }

        image.SaveAsBmp(path);
    }
    
    public static void WriteMatrixToFile(int[,] matrix, string path)
    {
        int rows = matrix.GetLength(0);
        int cols = matrix.GetLength(1);

        string[] lines = new string[rows];

        for (int i = 0; i < rows; i++)
        {
            string[] rowValues = new string[cols];
            for (int j = 0; j < cols; j++)
            {
                rowValues[j] = matrix[i, j].ToString();
            }

            lines[i] = string.Join(" ", rowValues);
        }

        File.WriteAllLines(path, lines);
    }


    public static bool CompareMatrix(int[,]? a, int[,]? b)
    {
        if (a == null || b == null) return false;
        int ah = a.GetLength(0);
        int aw = a.GetLength(1);
        int bh = b.GetLength(0);
        int bw = b.GetLength(1);
        if (ah != bh || aw != bw) return false;

        for (int y = 0; y < ah; y++)
        {
            for (int x = 0; x < aw; x++)
            {
                if (a[y, x] != b[y, x]) return false;
            }
        }

        return true;
    }

    public static int[,] ReadMatrixFromFile(string path)
    {
        string[] lines = File.ReadAllLines(path);
        var firstLine = lines[0].ToCharArray();

        int rows = lines.Length;
        int cols = firstLine.Length;

        int[,] matrix = new int[rows, cols];

        for (int i = 0; i < rows; i++)
        {
            var values = lines[i].ToCharArray();
            if (values.Length == firstLine.Length)
            {
                for (int j = 0; j < cols; j++)
                {
                    matrix[i, j] = values[j] == '1' ? 1 : 0;
                }
            }
            else
            {
                Console.WriteLine($"Line out of range {i}");
            }
        }

        return matrix;
    }
    
    public static void PrintMatrix(int[,] matrix)
    {
        int rows = matrix.GetLength(0);
        int cols = matrix.GetLength(1);

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                Console.Write(matrix[i, j] + " ");
            }

            Console.WriteLine();
        }
    }

    public static List<(int[,] signal, int errors)> DetectSignal(int[,] matrix)
    {
        var result = new List<(int[,] signal, int errors)>();
        int lastAcceptedX = int.MinValue;
        for (int x = 0; x < matrix.GetLength(1); x++)
        {
            int[,] test = new int[matrix.GetLength(0), matrix.GetLength(1)];
            bool detected = true;
            int errors = 0;

            for (int y = 0; y < matrix.GetLength(0); y++)
            {
                int sum = (x <= 0 ? 0 : matrix[y, x - 1])
                          + matrix[y, x + 0]
                          + (x >= matrix.GetLength(1) - 1 ? 0 : matrix[y, x + 1]);
                if (sum > 0)
                {
                    if (matrix[y, x - 0] == 0)
                    {
                        errors++;
                    }

                    if (x - 1 >= 0)
                    {
                        test[y, x - 1] = matrix[y, x - 1];
                    }

                    test[y, x + 0] = matrix[y, x + 0];
                    if (x < matrix.GetLength(1) - 1)
                    {
                        test[y, x + 1] = matrix[y, x + 1];
                    }
                }
                else
                {
                    detected = false;
                    break;
                }
            }

            if (detected && errors <= 0.5 * matrix.GetLength(0))
            {
            
                if (result.Count > 0 && x == lastAcceptedX + 1)
                {
                    if (errors < result[^1].errors)
                    {
                        result[^1] = (test, errors); 
                        lastAcceptedX = x;
                    }
                }
                else
                {
                    result.Add((test, errors));
                    lastAcceptedX = x;
                }
            }
        }

        return result;
    }
}