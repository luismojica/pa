import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

public class GroupGeneration {
	public static void main(String [] args) {
		String cmd = args[0];
		if(cmd.equalsIgnoreCase("agl"))
			print(AGL(Integer.parseInt(args[1]), Integer.parseInt(args[2])));
		else if(cmd.equalsIgnoreCase("pgl"))
			print(PGL(Integer.parseInt(args[1]), Integer.parseInt(args[2])));
		else if(cmd.equalsIgnoreCase("gfmult"))
			print(GF(Integer.parseInt(args[1]), Integer.parseInt(args[2]))[1]);
		else if(cmd.equalsIgnoreCase("gfadd"))
			print(GF(Integer.parseInt(args[1]), Integer.parseInt(args[2]))[0]);
		else if(cmd.equalsIgnoreCase("pgaml"))
			print(PGamL(Integer.parseInt(args[1]), Integer.parseInt(args[2])));
		else if(cmd.equalsIgnoreCase("agaml"))
			print(AGamL(Integer.parseInt(args[1]), Integer.parseInt(args[2])));
		else if(cmd.equalsIgnoreCase("aglcos"))
			printCosets(AGLCosets(Integer.parseInt(args[1]), Integer.parseInt(args[2])));
		else
			System.out.println("Invalid argument. Possible arguments: agl, pgl, gfmult, gfadd, pgaml, agaml, aglcos");
	}
	
	public static void printCosets(int[][][] blocks) {
		for(int i = 0; i < blocks.length; i++) {
			System.out.printf("\nCoset %d:\n", i+1);
			print(blocks[i]);
		}
	}
	
	public static int[][][] AGLCosets(int p, int r) {
		int[][] agl = AGL(p,r);
		int n = (int)Math.pow(p, r);
		int[][][] res = new int[n-1][n][];
		int block = 0;
		int row = 0;
		for(int perm = 0; perm < agl.length; perm++) {
			res[block][row++] = agl[perm];
			if(row == n) {
				block++;
				row = 0;
			}
		}
		return res;
	}
	
	public static void print(int[][] arr) {
		for(int i = 0; i < arr.length; i++) {
			for(int j = 0; j < arr[i].length; j++) {
				System.out.printf("%4d", arr[i][j]);
			}
			System.out.println();
		}
	}
	
	public static int[][] AGamL(int p, int r) {
		int n = (int)Math.pow(p, r);
		int[][] parr = new int[r*n*(n-1)][n];
		int[][] multTable = GF(p,r)[1];
		int k = 0;
		int[][] agl = AGL(p,r);
		for(int i = 0; i < r; i++) {
			int[] pwrPerm = powerMap(p, i, n, multTable);
			for(int perm = 0; perm < agl.length; perm++)
				parr[k++] = compose(agl[perm], pwrPerm);
		}
		return parr;
	}
	
	public static int[][] PGamL(int p, int r) {
		int n = (int)Math.pow(p, r);
		int[][] parr = new int[(n+1)*n*(n-1)*r][n+1];
		int[][] multTable = GF(p,r)[1];
		int k = 0;
		int[][] pgl = PGL(p,r);
		for(int i = 0; i < r; i++) {
			int[] pwrPerm = powerMap(p, i, n, multTable);
			for(int perm = 0; perm < pgl.length; perm++)
				parr[k++] = compose(pgl[perm], pwrPerm);
		}
		return parr;
	}
	
	public static int[] compose(int[] a, int[] b) {
		int[] res = new int[Math.max(a.length, b.length)];
		for(int i = 0; i < res.length; i++)
			if(i < b.length && i < a.length)
				res[i] = a[b[i]];
			else if(i < a.length)
				res[i] = a[i];
			else if(i < b.length)
				res[i] = b[i];
		return res;
	}
	
	public static int[] powerMap(int p, int i, int mod, int[][] mult) {
		int[] res = new int[mod];
		for(int n = 0; n < res.length; n++)
			res[n] = pow(n, (int)Math.pow(p, i), mult);
		return res;
	}
	
	public static int pow(int a, int r, int[][] mult) {
		if(r == 0)
			return 1;
		else if(r == 1)
			return a;
		else if(r % 2 == 0)
			return pow(mult[a][a], r/2, mult);
		else
			return mult[a][pow(mult[a][a], (r-1)/2, mult)];
	}
	
	public static int[][] PGL(int p, int r) {
		int n = (int)Math.pow(p, r);
		int[][] pArr = new int[(n+1)*n*(n-1)][n+1];
		int[][][] tables = GF(p,r);
		int k = 0;
		for(int a = 0; a < n; a++)
			for(int b = 0; b < n; b++)
				for(int c = 0; c < n; c++)
					for(int d = 0; d < n; d++) {
						if(tables[1][a][d] != tables[1][c][b] && k != pArr.length) {
							int[] temp = evaluatePGL(tables, a, b, c, d);
							if(!contains(pArr, temp, 0, k))
								pArr[k++] = temp;
						}
					}
		return pArr;
	}
	
	public static boolean contains(int[][] arr, int[] p, int s, int t) {
		for(int i = s; i < t; i++) {
			boolean same = true;
			for(int j = 0; j < p.length; j++)
				if(arr[i][j] != p[j]) {
					same = false;
					break;
				}
			if(same)
				return true;
		}
		return false;
	}
	
	public static int[] evaluatePGL(int[][][] tables, int a, int b, int c, int d) {
		int[] res = new int[tables[0].length+1];
		for(int i = 0; i < res.length-1; i++) {
			int den = tables[2][0][tables[0][tables[1][c][i]][d]];
			if(den == -1)
				res[i] = res.length-1;
			else
				res[i] = tables[1][den][tables[0][tables[1][a][i]][b]];
		}
		if(c == 0)
			res[res.length-1] = res.length-1;
		else
			res[res.length-1] = tables[1][a][tables[2][0][c]];
		return res;
	}
	
	public static int[][] AGL(int p, int r) {
		int n = (int)Math.pow(p, r);
		int[][] pArr = new int[n*(n-1)][n];
		int[][][] tables = GF(p,r);
		int k = 0;
		for(int i = 1; i < n; i++)
			for(int j = 0; j < n; j++)
				pArr[k++] = evaluateAGL(tables, i, j);
		return pArr;
	}
	
	public static int[] evaluateAGL(int[][][] tables, int a, int b) {
		int[] res = new int[tables[0].length];
		for(int i = 0; i < res.length; i++)
			res[i] = tables[0][tables[1][a][i]][b];
		return res;
	}
	
	public static int[][][] GF(int prime, int power) {
		if(!(new BigInteger(Integer.toString(prime))).isProbablePrime(10000)) {
			System.out.println("Not a prime. Please don't try to break me.");
			System.exit(1);
		}
		Polynomial.mod = prime;
		Polynomial irr = findRandomIrreducible(prime, power);
		int[][][] res = new int[3][][];
		ArrayList<Polynomial> arr = genPolynomials(prime, power);
		
		res[0] = addTable(prime, power, arr);
		res[1] = multTable(prime, power, irr, arr);
		res[2] = new int[1][];
		res[2][0] = inverseTable(res[1], arr);
		return res;
	}
	
	public static int[][] multTable(int prime, int pow, Polynomial irr, ArrayList<Polynomial> arr) {
		Polynomial.mod = prime;
		int sz = (int)Math.pow(prime, pow);
		int[][] table = new int[sz][sz];
		for(int a = 0; a < sz; a++)
			for(int b = 0; b < sz; b++) {
				Polynomial result = arr.get(a).mult(arr.get(b)).divide(irr)[1];
				for(int i = 0; i < arr.size(); i++)
					if(arr.get(i).equals(result)) {
						table[a][b] = i;
						break;
					}
			}
		return table;
	}
	
	public static int[][] addTable(int prime, int pow, ArrayList<Polynomial> arr) {
		Polynomial.mod = prime;
		int sz = (int) Math.pow(prime, pow);
		int[][] table = new int[sz][sz];
		for(int a = 0; a < sz; a++)
			for(int b = 0; b < sz; b++) {
				Polynomial result = arr.get(a).add(arr.get(b));
				for(int i = 0; i < arr.size(); i++) {
					if(arr.get(i).equals(result)) {
						table[a][b] = i;
						break;
					}
				}
			}
		return table;
	}
	
	public static int[] inverseTable(int[][] mult, ArrayList<Polynomial> arr) {
		int one = -1;
		int[] inv = new int[arr.size()];
		for(int i = 0; i < arr.size(); i++)
			if(arr.get(i).deg == 0 && arr.get(i).coef[0] == 1)
				one = i;
		for(int i = 1; i < mult.length; i++)
			for(int j = 0; j < mult[i].length; j++)
				if(mult[i][j] == one)
					inv[i] = j;
		inv[0] = -1;
		return inv;
	}
	
	public static ArrayList<Polynomial> genPolynomials(int mod, int deg) {
		ArrayList<Polynomial> res = new ArrayList<Polynomial>();
		int[] coef = new int[deg];
		recurse(coef, 0, res, mod);
//		for(int i = 1; i < res.size(); i++)
//			if(res.get(i).deg == 0 && res.get(i).coef[0] == 1) {
//				Polynomial temp = res.get(1);
//				res.set(1, res.get(i));
//				res.set(i, temp);
//				break;
//			}
		return res;
	}
	
	public static int[] reverse(int[] arr) {
		int[] res = new int[arr.length];
		for(int i = 0; i < res.length; i++)
			res[i] = arr[arr.length-1-i];
		return res;
	}
	
	public static void recurse(int[] coef, int i, ArrayList<Polynomial> res, int mod) {
		if(i == coef.length) {
			res.add(new Polynomial(reverse(coef)));
			return;
		}
		for(int k = 0; k < mod; k++) {
			coef[i] = k;
			recurse(coef, i+1, res, mod);
		}
	}
	
	public static Polynomial findRandomIrreducible(int prime, int power) {
		int[] temp = new int[(int)Math.pow(prime, power) + 1];
		temp[temp.length-1] = 1;
		temp[1] = prime - 1;
		Polynomial cur = randomPolynomial(prime, power);
		while(isReducible(cur, prime))
			cur = randomPolynomial(prime, power);
		return cur;
	}
	
	public static boolean isReducible(Polynomial p, int mod) {
		for(int i = 1; i < p.deg; i++) {
			int[] coef = new int[(int)Math.pow(mod, i)+1];
			coef[coef.length-1] = 1;
			coef[1] = mod-1;
			Polynomial test = new Polynomial(coef);
			if(Polynomial.gcd(p, test).deg > 0) {
				return true;
			}
		}
		return false;
	}
	
	public static Polynomial randomPolynomial(int mod, int deg) {
		int[] coef = new int[deg+1];
		Random rand = new Random();
		coef[coef.length-1] = 1;
		for(int i = 0; i < coef.length-1; i++)
			coef[i] = rand.nextInt(mod);
		return new Polynomial(coef);
	}
}

class Polynomial {
	static int mod = 0;
	int[] coef;
	int deg;
	
	public Polynomial(int[] coef) {
		deg = coef.length-1;
		while(deg >= 0 && coef[deg] == 0)
			deg--;
		this.coef = Arrays.copyOf(coef, deg+1);
	}
	
	public Polynomial mult(Polynomial other) {
		if(deg == -1 || other.deg == -1)
			return new Polynomial(new int[]{});
		int[] newCoef = new int[deg + other.deg + 1];
		for(int i = 0; i < coef.length; i++)
			for(int j = 0; j < other.coef.length; j++)
				newCoef[i+j] = (coef[i]*other.coef[j] + newCoef[i+j]) % mod;
		return new Polynomial(newCoef);
	}
	
	public Polynomial add(Polynomial other) {
		int[] newCoef = new int[Math.max(coef.length, other.coef.length)];
		for(int i = 0; i < newCoef.length; i++) {
			if(i < coef.length)
				newCoef[i] = (newCoef[i] + coef[i]) % mod;
			if(i < other.coef.length)
				newCoef[i] = (newCoef[i] + other.coef[i]) % mod;
		}
		return new Polynomial(newCoef);
	}
	
	public Polynomial subtract(Polynomial other) {
		int[] newCoef = new int[Math.max(coef.length, other.coef.length)];
		for(int i = 0; i < newCoef.length; i++) {
			if(i < coef.length)
				newCoef[i] = (coef[i]+ newCoef[i]) % mod;
			if(i < other.coef.length)
				newCoef[i] = (newCoef[i] - other.coef[i] + mod) % mod;
		}
		return new Polynomial(newCoef);
	}
	
	public Polynomial divideLeadTerms(Polynomial other) {
		int temp = deg - other.deg;
		int value = coef[coef.length-1] * invert(other.coef[other.coef.length-1]);
		int newCoef[] = new int[temp+1];
		newCoef[temp] = value;
		return new Polynomial(newCoef);
	}
	
	public Polynomial[] divide(Polynomial divisor) {
		Polynomial[] res = new Polynomial[2];
		if(divisor.deg == -1)
			return null;
		res[0] = new Polynomial(new int[]{});
		res[1] = this.copy();
		while(res[1].deg != -1 && res[1].deg >= divisor.deg) {
			Polynomial temp = res[1].divideLeadTerms(divisor);
			res[0] = res[0].add(temp);
			res[1] = res[1].subtract(temp.mult(divisor));
		}
		return res;
	}
	
	public static int invert(int n) {
		for(int a = 1; a < mod; a++)
			if((a * n) % mod == 1)
				return a;
		return -1;
	}
	
	public static Polynomial gcd(Polynomial a, Polynomial b) {
		if(b.deg == -1)
			return a;
		else
			return gcd(b, a.divide(b)[1]);
	}
	
	public Polynomial copy() {
		int[] newCoef = Arrays.copyOf(coef, coef.length);
		return new Polynomial(newCoef);
	}
	
	public boolean equals(Object other) {
		if(other instanceof Polynomial) {
			Polynomial poly = (Polynomial)other;
			if(poly.deg != deg)
				return false;
			for(int i = 0; i < coef.length; i++)
				if(coef[i] != poly.coef[i])
					return false;
			return true;
		}
		return false;
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		if(deg == -1)
			return "0";
		for(int i = coef.length-1; i >= 2; i--) {
			if(coef[i] == 1)
				sb.append("x^" + i + " + ");
			else if(coef[i] != 0)
				sb.append(coef[i] + "x^" + i + " + ");
		}
		if(coef.length >= 2 && coef[1] != 0)
			sb.append((coef[1] == 1 ? "" : coef[1]) + "x + ");
		if(coef.length >= 1 && coef[0] != 0)
			sb.append(coef[0]);
		else
			sb.delete(sb.length()-3, sb.length());
		return sb.toString();
	}
}