module ode.plot;

import std.traits: isIntegral;

enum DataFileFormat{
	csv,
	binaryMatrix,
	binaryGeneral
}

/***
 * Data container
 *****/
struct ScatterPlot(R, int Order)
if(Order == 1 || Order == 2){
	this(in string label) @safe pure nothrow{
		this._header= label;
	}

	// mutable
	@safe pure nothrow{
		static if(Order == 1){
			void append(in R t, in R x){
				this._t ~= t;
				this._x ~= x;
			}
		}
		else{
			void append(in R t, in R x, in R v){
				this._t ~= t;
				this._x ~= x;
				this._v ~= v;
			}
		}

		void setHeader(in string label){
			_header= label;
		}
	}

	/**********
	 * 
	 **/
	string output(DataFileFormat FileType)() @system{
		import std.array: appender;
		import std.format: formattedWrite;
		import std.range: lockstep;

		static if(FileType is DataFileFormat.csv){
			/***
			 * CSV
			 ***/
			auto buf= appender!string;
			if(label.length > 0) buf.writeln(_header);

			static if(Order == 1){
				foreach(R t, R x; lockstep(_t, _x)){
					buf.formattedWrite!"%f,%f\n"(t, x);
				}
			}
			else{
				foreach(R t, R x, R v; lockstep(_t, _x, _v)){
					buf.formattedWrite!"%f,%f,%f\n"(t, x, v);
				}
			}
			return buf.data;
		}
		else if(FileType is DataFileFormat.binaryMatrix || 
		FileType is DataFileFormat.binaryGeneral){
			/***
			 * binary matrix
			 ***/
			char[] buf;
			static if(Order == 1){
				foreach(R t, R x; lockstep(_t, _x)){
					this.num[]= [cast(float)t, cast(float)x];
					buf ~= this.data;
				}
			}
			else{
				foreach(R t, R x, R v; lockstep(_t, _x, _v)){
					this.num[]= [cast(float)t, cast(float)x, cast(float)v];
					buf ~= this.data;
				}
			}
			return buf.idup;
		}
/+		else if(FileType is DataFileFormat.binaryGeneral){
			/***
			 * binary general
			 ***/
		}+/
		else{}

		//return buf.data;
	}


protected:
	string _header;
	static if(Order == 1){
		union{
			float[2] num;
			char[8] data;
		};
	}
	else{
		union{
			float[3] num;
			char[12] data;
		};
	}

	R[] _t;
	R[] _x;
	static if(Order == 2){
		R[] _v;
	}
}
