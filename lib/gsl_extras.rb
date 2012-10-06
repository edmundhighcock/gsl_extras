####################################
# = GSL Tools                    
####################################
#
# A set of tools which extend the Ruby wrapper for the GNU Scientific Library



#######################
# Fixes for rb-gsl    #
#######################

# Means that #inspect, eval, Marshal.dump etc all work in the expected way

	class GSL::Vector
# 		aliold :join
		def to_gslv
			self
		end
		def join(*args)
			to_a.join(*args)
		end
		def inspect
			"GSL::Vector.alloc(#{to_a.inspect})"
		end
	end
	class GSL::Vector
		def _dump(depth)
			return Marshal.dump(self.to_a)
		end
		def self._load(string)
			return self.alloc(Marshal.load(string))
		end
	end
	class GSL::Matrix
		def _dump(depth)
			return Marshal.dump(self.to_a)
		end
		def self._load(string)
			arr = Marshal.load(string)
			return self.alloc(arr.flatten, arr.size, arr[0].size)
		end
		def inspect
			arr=self.to_a
			"GSL::Matrix.alloc(#{arr.flatten.inspect}, #{arr.size}, #{arr[0].size})"
		end
	end

	class GSL::Complex
		def to_a
			return [self.real, self.imag]
		end
		def _dump(depth)
			return Marshal.dump(self.to_a)
		end
		def self._load(string)
			return self.alloc(Marshal.load(string))
		end
		def inspect
			"GSL::Complex.alloc(#{self.to_a.inspect})"
		end
	end
	class GSL::Vector::Complex
		def _dump(depth)
			return Marshal.dump([self.real.to_a, self.imag.to_a])
		end
		def self._load(string)
			re, im = Marshal.load(string)
			return self.alloc(re.zip(im))
		end
		def inspect
			re, im = [self.real.to_a, self.imag.to_a]
			"GSL::Vector::Complex.alloc(#{re.zip(im).inspect})"
		end

	end		



class CodeRunner
	NaN = GSL::NAN
	Infinity = GSL::POSINF
end



module GSL

# A class for interpolation from scattered multidimensional data using radial basis functions.
#
# E.g. 
# 	x = GSL::Vector.alloc([0,3,6])
# 	y = GSL::Vector.alloc([0,1,2])
# 	z = GSL::Vector.alloc([0,5,10])
# 	normalising_radius = GSL::Vector.alloc([3,1])
#
# 	int = GSL::ScatterInterp.alloc(:linear, [x,y,z], false, normalising_radius)
# 	puts int.eval(4.5, 1.7)

	
	
class ScatterInterp
	
	THIN_PLATE_SPLINES = :thin_plate_splines
	
	# Create a new interpolation class, for interpolating a function of one or more variables from a scattered dataset. datavecs is an array of vectors; the last vector should be the values of the function to be interpolated from, the other vectors are the coordinates or parameters corresponding to those values. Norm is a boolean; should the normalised functions be used? Default false. Func is the particular basis function to be used: can be
	# * :linear
	# * :cubic
	# * :thin_plate_splines
	# * :multiquadratic
	# * :inverse_multiquadratic
	#
	# <tt>r0</tt> is a normalising radius which should be on the order of the scale length of the variation. If the scale length differs for each direction, an array of values can be passed; the length of this array should be one less than the number of datavecs, i.e. it should be the number of coordinates.
	
	def self.alloc(func, datavecs, norm, r0=1.0)
		new(func, datavecs, norm, r0)
	end
	
	
	
	def initialize(func, datavecs, norm=false, r0=1.0)
# 		p datavecs
		 @norm = norm
		@npoints = datavecs[0].size
		datavecs.map!{|v| v.to_a}
		datavecs.each{|vec| raise ArgumentError.new("Datavectors must all have the same size ") unless vec.size == @npoints}
		@func = func
		data = datavecs.pop
		@gridpoints = GSL::Matrix.alloc(*datavecs).transpose
# 			puts @gridpoints.shape
		@dim = datavecs.size
		@r0 = r0
		if @r0.kind_of? Numeric
			v = GSL::Vector.alloc(@dim)
			v.set_all(@r0)
			@r0 = v
		end
		m = GSL::Matrix.alloc(@npoints, @npoints)
		for i in 0...@npoints
			for j in 0...@npoints
# 					ep i, j
# 					if true or i>= j
# 					p @gridpoints.row(i), @gridpoints.row(j) if i == j
				m[i,j] = function(@gridpoints.row(i), @gridpoints.row(j)) #(radius(@gridpoints.row(i), @gridpoints.row(j)))
# 					else 
# 						m[i,j] = 0.0
# 					end
			end
# 				data[i] = data[i] * m.get_row(i).sum if norm
		end
# 			ep m
		@weights = m.LU_solve(GSL::Vector.alloc(data))
# 			ep @weights
	end
	
	
	def radius(vec1, vec2) 	# :nodoc: 
		Math.sqrt((vec1 -  vec2).square.sum)
	end
	
	# :nodoc:
	
	def normalized_radius(vec1, vec2)
		#Math.sqrt(((vec1 - vec2) / @r0).square.sum)
		case @r0
		when Numeric
			Math.sqrt(((vec1 - vec2) / @r0).square.sum)
		else 
			Math.sqrt(((vec1/@r0.to_gslv - vec2/@r0.to_gslv)).square.sum)
		end

	end
	
	# Return the value of the interpolation kernel for the separation between the two given vectors. If linear was chosen this will just be the normalised distance between the two points.
	
	def function(vec1, vec2)
		case @func
		when :linear
			return normalized_radius(vec1, vec2)
		when :cubic_alt
			return normalized_radius(vec1, vec2)**(1.5)
		when :thin_plate_splines
			return 0.0 if radius(vec1, vec2) == 0.0
			return normalized_radius(vec1, vec2)**2.0 * Math.log(normalized_radius(vec1, vec2))
		when :thin_plate_splines_alt
			rnorm = ((@r0.prod.abs)**(2.0/@r0.size)*(((vec1-vec2).square / @r0.square).sum + 1.0))
			return rnorm * Math.log(rnorm)
		when :multiquadratic
# 			return Math.sqrt(radius(vec1, vec2)**2 + Math.sqrt(@r0.square.sum))
			(@r0.prod.abs)**(1.0/@r0.size)*Math.sqrt(((vec1-vec2).square / @r0.square).sum + 1.0)
		when :inverse_multiquadratic
			1.0 / ((@r0.prod.abs)**(1.0/@r0.size)*Math.sqrt(((vec1-vec2).square / @r0.square).sum + 1.0))
		when :cubic
			((@r0.prod.abs)**(2.0/@r0.size)*(((vec1-vec2).square / @r0.square).sum + 1.0))**(1.5)
# 			invs = ((vec1-vec2).square + @r0.square).sqrt**(-1)
# 			invs.sum
# 				p @ro
# 			return 1.0 / Math.sqrt(radius(vec1, vec2)**2 + Math.sqrt(@r0.square.sum))
# 		when :inverse_multiquadratic
# # 				p @ro
# 			return 1.0 / Math.sqrt(radius(vec1, vec2)**2 + Math.sqrt(@r0.square.sum))
		else
			raise ArgumentError.new("Bad radial basis function: #{@func}")
		end
	end
	
	# Return the interpolated value for the given parameters.
	
	def eval(*pars)
		raise ArgumentError("wrong number of points") if pars.size != @dim
# 			p vals
		pars = GSL::Vector.alloc(pars)
		return @npoints.times.inject(0.0) do |sum, i|
# 			sum + function(radius(vals, @gridpoints.row(i)))*@weights[i]
			sum + function(pars, @gridpoints.row(i))*@weights[i]
		end
	end
	
	# Evaluate the function, 
	
	def gaussian_smooth_eval(*vals, sigma_vec)
		npix = 7
		raise "npix must be odd" if npix%2==0
		case vals.size
		when 2
# 			delt0 = 3.0*0.999999*sigma_vec[0] / (npix-1)
# 			delt1 = 3.0*0.999999*sigma_vec[1] / (npix-1)
# 			sig3 = 3.0*sigma
			vals0 = GSL::Vector.linspace(vals[0] - 3.0* sigma_vec[0], vals[0] + 3.0* sigma_vec[0], npix)
			vals1 = GSL::Vector.linspace(vals[1] - 3.0* sigma_vec[1], vals[1] + 3.0* sigma_vec[1], npix)
			mat = GSL::Matrix.alloc(vals0.size, vals1.size)
			for i in 0...vals0.size
				for j in 0...vals1.size
					mat[i,j] = eval(vals0[i], vals1[j])
				end
			end
			mat.gaussian_smooth(*sigma_vec.to_a)
			cent = (npix - 1) / 2
			return mat[cent, cent]
			
		else
			raise 'Not supported for this number of dimensions yet'
		end
	end
		
	# Create a GSL::Contour object for making contours of the interpolated function. Only works for functions of 2 variables.
		
	def to_contour(grid_size=@npoints)
		m = Matrix.alloc(grid_size, grid_size)
		raise TypeError("Must be 3d data") unless @gridpoints.shape[1] == 2
# 			p @gridpoints.shape
# 			p @gridpoints
		xmax, xmin = @gridpoints.col(0).max, @gridpoints.col(0).min
		ymax, ymin = @gridpoints.col(1).max, @gridpoints.col(1).min
		p 'x', xmax, 'y', ymax
		xvec = GSL::Vector.alloc((0...grid_size).to_a) * (xmax - xmin) / (grid_size - 1.0).to_f + xmin
		yvec = GSL::Vector.alloc((0...grid_size).to_a) * (ymax - ymin) / (grid_size - 1.0).to_f + ymin	
		p 'x', xvec.max, 'y', yvec.max

		for i in 0...grid_size
			for j in 0...grid_size
# 					p xvec[i], yvec[j]
				m[i,j] = eval(xvec[i], yvec[j])
			end
		end
		p 'm', m.max
		Contour.alloc(xvec, yvec, m)
	end
	
end


# A class for making contours of a function on a regular two dimensional grid. If contours of scattered data are required, see GSL::ScatterInterp#to_contour.

class Contour
	
	# Create a new Contour object. <tt>x</tt> and <tt>y</tt> are vectors of coordinates, and <tt>grid</tt> is a matrix of values on those coordinates.
	
	def self.alloc(x, y, grid)
		new(x, y, grid)
	end
	
	attr_accessor :keep_path_data

	def initialize(x, y, grid)
		@x = x; @y=y; @grid=grid
# 			p @grid, @x, @y
		raise ArgumentError.new("Unmatching data sizes: #{x.size}, #{y.size}, #{grid.shape}") unless [x.size, y.size] == grid.shape
		@adaptive = false			
	end
	
	
	
	def set_adaptive(func, scale, multi_adaptive=false) # :nodoc:
		@func = func; @adaptive = true
		@multi_adaptive = multi_adaptive
		if @multi_adaptive
			@adaption_scale = 4
			raise "Adaption scale should be a power of two for multi_adaptive contour generation" if scale % 2 == 1 and not scale == 1
		
			@next_adaption_scale = scale / 2
		else
			@adaption_scale = scale*2
		end
	end
	
	# Create a series of contours at the given values. Returns a hash of {value => array_of_contours}. The array_of_contours is a an array of arrays, where each array is a list of [x, y] coordinates along the contour.
	
	def contours(*values)
		(values = (0..values[0]+1).to_a.map{|i| i.to_f * (@grid.max - @grid.min) / ( values[0]+1) + @grid.min}; values.pop; values.shift) if values.size==1 and values[0].kind_of? Integer
		cons = values.inject({}){|hash, val| hash[val] = []; hash}
# 			p cons
		for i in 0...((@x.size / 2.0).ceil - 1)
			for j in 0...((@y.size / 2.0).ceil - 1)
				analyse_cell(i*2, j*2, cons)
			end
		end
# 			pp cons
		cons.keys.each{|val| cons[val] = connect_contours(cons[val])}
		@last_contours = cons
	end
	
	# Create a GraphKit object of the contours.
	
	def graphkit(*args)
		if args.size == 0
			conts = @last_contours
		else
			conts = contours(*args)
		end
		graphs = conts.map do |val, cons|
			unless cons[0]
				nil
			else
				(cons.map do |con|
	# 				p con
					contour = con.transpose
					kit = CodeRunner::GraphKit.autocreate({x: {data: contour[0]}, y: {data: contour[1], title: val.to_s}})
					kit.data[0].with = "l"
					kit
				end).sum
			end
		end
		graphs.compact.reverse.sum
	end
		
	
	#Edges:       1 __
	#         0 | 4 \ /5  | 2
	#           | 7 / \ 6 |
	#             3 __
	
	VALID_CONNECTIONS = [[0,4,5,2], [0,4,1], [0,7,3], [0,7,6,2], [0,4,5,6,3],[0,7,6,5,1], [1,4,7,3], [1,5,6,3], [1,5,2], [1,4,7,6,2], [2,6,3],[2,5,4,7,3]]
	
	def get_crossed_edges(i, j, cons)
		ce = {}
		edges = {0=>[[i, j], [i+2, j]],
				1=>[[i, j], [i, j+2]],
				2=>[[i, j+2], [i+2, j+2]],
				3=>[[i+2, j], [i+2, j+2]],
				4=>[[i, j], [i+1, j+1]],
				5=>[[i, j+2], [i+1, j+1]],
				6=>[[i+2, j+2], [i+1, j+1]],
				7=>[[i+2, j], [i+1, j+1]]
			}    
		cons.keys.each do |val|
			ce[val] = {}
			edges.each do |edge, (start, fin)|
				bounds = [@grid[*start], @grid[*fin]]
# 					p edge, bounds if bounds.max < 4
				if val <= bounds.max and val > bounds.min
					dx = @x[fin[0]] - @x[start[0]]
					dy = @y[fin[1]] - @y[start[1]]
					df = bounds[1] - bounds[0]
					xcross = @x[start[0]] + (val-bounds[0]) / df * dx
					ycross = @y[start[1]] + (val-bounds[0]) / df * dy
					ce[val][edge] = [xcross, ycross]
				end
				
			end
		end
		ce
	end
		
	private :get_crossed_edges

	def analyse_cell(i, j, cons)
		crossed_edges = get_crossed_edges(i, j, cons)
		if @adaptive and crossed_edges.values.find_all{|crossings| crossings.size>0}.size > 0
			@adaptive_matrix ||= GSL::Matrix.alloc(@adaption_scale+1, @adaption_scale+1)
			@adaptive_x ||= GSL::Vector.alloc(@adaption_scale+1)
			@adaptive_y ||= GSL::Vector.alloc(@adaption_scale+1)
			dx = (@x[i+2] - @x[i])/@adaption_scale.to_f
			dy = (@y[j+2] - @y[j])/@adaption_scale.to_f
			x = @x[i]
			y = @y[j]
			for ii in 0...@adaption_scale+1
				for jj in 0...@adaption_scale+1
					@adaptive_x[ii] = x + ii * dx
					@adaptive_y[jj] = y + jj * dy
					@adaptive_matrix[ii, jj] = @func.eval(@adaptive_x[ii], @adaptive_y[jj])
				end
			end
			cell_contour = Contour.alloc(@adaptive_x, @adaptive_y, @adaptive_matrix)
			cell_contour.keep_path_data = true
			if @multi_adaptive and @next_adaption_scale > 1
				#p @next_adaption_scale
				#p "#{@adaptive_x.max}, #{@adaptive_x.min}, #{@x[i]}, #{@x[i+2]}"
				cell_contour.set_adaptive(@func, @next_adaption_scale, true)
			end

			cell_cons = cell_contour.contours(*cons.keys)
			#kit = cell_contour.graphkit
			#kit.gnuplot(xrange: [@x[i], @x[i+2]], yrange: [@y[j], @y[j+2]], style: "data lp") 
			#p cell_cons
			cell_cons.each do |val, cell_val_contours|
				cell_val_contours.each do |path|
					# we have to relabel the path from the fine scale contour so that it fits with the course scale labelling system. only the beginning and end of the path matters
					path_start = path[0][0]  
					if path_start[0] == @x[i]
							path[0][1] = [i, j, 1]
					elsif path_start[0] == @x[i+2]
						  path[0][1] = [i, j, 3]
					else
						if path_start[1] == @y[j]
						  path[0][1] = [i, j, 0]
						elsif path_start[1] == @y[j+2]
						  path[0][1] = [i, j, 2]
						else
							raise "Could not find path_start; #{path_start}; x: #{@x[i]}, #{@x[i+2]}; y: #{@y[j]}, #{@y[j+2]} "
						end
					end
					path_end = path[-1][0]  
					if path_end[0] == @x[i]
							path[-1][1] = [i, j, 1]
					elsif path_end[0] == @x[i+2]
						  path[-1][1] = [i, j, 3]
					else
						if path_end[1] == @y[j]
						  path[-1][1] = [i, j, 0]
						elsif path_end[1] == @y[j+2]
						  path[-1][1] = [i, j, 2]
						else
							raise "Could not find path_end #{path_end}; x: #{@x[i]}, #{@x[i+2]}; y: #{@y[j]}, #{@y[j+2]} "
						end
					end

					cons[val].push path
				end
			end
			#kit.close
		else
			crossed_edges.each do |val, crossings|
				outer = crossings.keys.find_all{|edge| edge < 4}
				inner = crossings.keys.find_all{|edge| edge > 4}
				next if outer.size == 0 and inner.size == 0
				VALID_CONNECTIONS.each do |connection|
					path = crossings.values_at(*connection).compact.zip(connection.map{|edge| [i, j, edge]})
	# 					p path
					next if path.size != connection.size
					connection.each{|edge| crossings.delete(edge)}
					cons[val].push path
					
				end
			end

# 				p val, crossings, cons
		end
	end
	
	private :analyse_cell
	
	def connect_contours(contours)
		return contours if contours.size == 1
		loop do
			catch(:restart) do 
# 				old_contours = contours
# 				contours = []
			joined = []
			for i in 0...contours.size
				break if i >= contours.size
				for j in i+1...contours.size
					break if j >= contours.size
					coni = contours[i]
					conj = contours[j]
					if joins?(coni[-1], conj[0])
						contours[i] = coni + conj
						contours.delete_at(j)
# 							redo
						throw(:restart)
					elsif joins?(coni[0], conj[-1])
						contours[i] = conj + coni
						contours.delete_at(j)
# 							redo
						throw(:restart)
					elsif joins?(coni[0], conj[0])
						contours[i] = coni.reverse + conj
						contours.delete_at(j)
# 							redo
						throw(:restart)
					elsif joins?(coni[-1], conj[-1])
						contours[i] = conj + coni.reverse
						contours.delete_at(j)
# 							redo
						throw(:restart)
					end


					
				end
			end
			contours.each{|con| con.uniq!}
			contours.map!{|con| con.map{|point| point[0]}} unless @keep_path_data
			return contours
			end
			
		end
	end
	
	private :connect_contours
			
	
	
	def joins?(path, contour)	
# 			p @x.size
# 			dx = (@x.subvector(1, @x.size - 1) - @x.subvector(0, @x.size - 1)).sum / @x.size #The average dx between gridpoints
# 			dy = (@y.subvector(1, @y.size - 1) - @y.subvector(0, @y.size - 1)).sum / @y.size
# 			eps = 1.0
# # 			return true
# 			return ((path[0] - contour[0]).abs < dx.abs * eps and (path[1] - contour[1]).abs < dy.abs * eps)
# 			p dx; exit
		pi, pj, pedge = path[1]
		ci, cj, cedge = contour[1]
		joins = false
		joins = ((pi == ci and pj == cj + 2 and pedge == 0 and cedge == 2) or
			(pi == ci and pj == cj - 2 and pedge == 2 and cedge == 0) or
			(pi == ci + 2 and pj == cj  and pedge == 1 and cedge == 3) or
			(pi == ci - 2 and pj == cj  and pedge == 3 and cedge == 1))
# 			unless joins
# 				p path[1], contour[1]
# 				STDIN.gets
# 			end
# 			if path[1] == [2,10,1]
# 				p path[1], contour[1], joins
# 				STDIN.gets
# 			end
		return joins
	end
	
	private :joins?
		
		

	
end
class Contour2
	
	# Create a new Contour object. <tt>x</tt> and <tt>y</tt> are vectors of coordinates, and <tt>grid</tt> is a matrix of values on those coordinates. If a function is given, it will be used to evaluate gridpoints rather than the matrix (though a matrix must still be provided). Provide a function if your contours only cover a small part of the domain and your function is expensive to evaluate.
	
	def self.alloc(x, y, grid, function=nil)
		new(x, y, grid)
	end
	
	
	def initialize(x, y, grid, function=nil)
		@function = function
		@evaluated = {}
		@x = x; @y=y; @grid=grid
# 			p @grid, @x, @y
		raise ArgumentError.new("Unmatching data sizes: #{x.size}, #{y.size}, #{grid.shape}") unless [x.size, y.size] == grid.shape
		@adaptive = false			
	end
	
	
	
	def set_adaptive(func) # :nodoc:
		@func = func; @adaptive = true
	end
	
	# Create a series of contours at the given values. Returns a hash of {value => array_of_contours}. The array_of_contours is a an array of arrays, where each array is a list of [x, y] coordinates along the contour.
	
	def contours(*values)
		(values = (0..values[0]+1).to_a.map{|i| i.to_f * (@grid.max - @grid.min) / ( values[0]+1) + @grid.min}; values.pop; values.shift) if values.size==1 and values[0].kind_of? Integer
		cons = values.inject({}){|hash, val| hash[val] = []; hash}
		get_startpoints(cons)
		#@analysed = {}
		@found = {}
		#p cons
		starts = cons

		cons = values.inject({}){|hash, val| hash[val] = []; hash}
		temp_cons = values.inject({}){|hash, val| hash[val] = []; hash}
		starts.each do |val, arr|
			#p 'arr', arr
			arr.each do |start_con|
			#arr.map{|edges| edges[-1][1].slice(0..1)}.each do |starti, startj|
				starti, startj = start_con[-1][1].slice(0..1)	
				temp_cons[val] = [start_con]
				p 'startj',  starti, startj
				loop do
					starti, startj = trace_contour(val, temp_cons, starti, startj)
					break unless starti and startj
				end
				cons[val].push temp_cons[val][0]
			end
		end
		cons.keys.each{|val| cons[val] = connect_contours(cons[val]);
			cons[val].map!{|con| con.map{|point| point[0]}}}
		
		@last_contours = cons
		gk =  graphkit
		#gk.gp.style = 'data with linespoints'
		gk.data.each{|dk| dk.with = 'lp'}
		gk.gnuplot
		@last_contours
	end

	def trace_contour(val, cons, starti, startj)
				old_start = cons[val][0][0][1]; old_end = cons[val][0][-1][1]
				#p 'old_start', old_start, 'old_end', old_end, 'starti', starti, 'startj', startj
				for delti in -1..1
					for deltj in -1..1
						celli, cellj = [starti + 2 * delti, startj + deltj * 2]
							#unless (@analysed[val] and @analysed[val][[celli, cellj]] ) or celli > (@x.size-3) or cellj >(@y.size - 3) or celli < 0 or cellj < 0
							unless  celli > (@x.size-3) or cellj >(@y.size - 3) or celli < 0 or cellj < 0

					 		#p 'analysing', celli, cellj 
					 		analyse_cell(celli, cellj, cons, val) 
							end

					end
				end
					cons[val] = connect_contours(cons[val])
				#p cons[val]
					new_contour =  cons[val].find{|cont| old_start ==  cont[0][1] or old_end == cont[-1][1]}
					unless new_contour
						cons[val].map!{|con| con.reverse}
					new_contour =  cons[val].find{|cont| old_start ==  cont[0][1] or old_end == cont[-1][1]}
					end
					raise "no new contour" unless new_contour
					cons[val] = [new_contour]
				#p cons[val]
				#newtails  = cons[val].map{|con| con[-1]}
				if cons[val][0][0][1] != old_start 
					return cons[val][0][0][1].slice(0..1)
				elsif cons[val][0][-1][1] != old_end
					return cons[val][0][-1][1].slice(0..1)
				else
					#p cons[val]
					
				p 'old_start', old_start, 'old_end', old_end
					#raise "Could not connect"
					return nil, nil
				end
	end

	def get_startpoints(cons)	
		[0,((@x.size / 2).floor - 1)-1].each do |i|
			for j in 0...((@y.size / 2).floor - 1)
				analyse_cell(i*2, j*2, cons)
			end
		end
		for i in 0...((@x.size / 2).floor - 1)
		[0,((@y.size / 2).floor - 1)-1].each do |j|
				analyse_cell(i*2, j*2, cons)
			end
		end
	end



	# Create a GraphKit object of the contours.
	
	def graphkit(*args)
		if args.size == 0
			conts = @last_contours
		else
			conts = contours(*args)
		end
		graphs = conts.map do |val, cons|
			unless cons[0]
				nil
			else
				(cons.map do |con|
	# 				p con
					contour = con.transpose
					kit = CodeRunner::GraphKit.autocreate({x: {data: contour[0]}, y: {data: contour[1], title: val.to_s}})
					kit.data[0].with = "l"
					kit
				end).sum
			end
		end
		graphs.compact.reverse.sum
	end
		
	
	#Edges:       1 __
	#         0 | 4 \ /5  | 2
	#           | 7 / \ 6 |
	#             3 __
	
	VALID_CONNECTIONS = [[0,4,5,2], [0,4,1], [0,7,3], [0,7,6,2], [0,4,5,6,3],[0,7,6,5,1], [1,4,7,3], [1,5,6,3], [1,5,2], [1,4,7,6,2], [2,6,3],[2,5,4,7,3]]
	
	def get_crossed_edges(i, j, cons, specific_val = nil)
		ce = {}
		edges = {0=>[[i, j], [i+2, j]],
				1=>[[i, j], [i, j+2]],
				2=>[[i, j+2], [i+2, j+2]],
				3=>[[i+2, j], [i+2, j+2]],
				4=>[[i, j], [i+1, j+1]],
				5=>[[i, j+2], [i+1, j+1]],
				6=>[[i+2, j+2], [i+1, j+1]],
				7=>[[i+2, j], [i+1, j+1]]
			}    
		cons.keys.each do |val|
			next if specific_val and specific_val != val
			ce[val] = {}
			edges.each do |edge, (start, fin)|
				if @function
					#p start, fin
					(@evaluated[start]= true, (@grid[*start] = @function.eval(@x[start[0]], @y[start[1]]))) unless @evaluated[start]
					(@evaluated[fin]= true, (@grid[*fin] = @function.eval(@x[fin[0]], @y[fin[1]])))  unless @evaluated[fin]
				end
				bounds = [@grid[*start], @grid[*fin]]
# 					p edge, bounds if bounds.max < 4
				if val <= bounds.max and val > bounds.min
					dx = @x[fin[0]] - @x[start[0]]
					dy = @y[fin[1]] - @y[start[1]]
					df = bounds[1] - bounds[0]
					xcross = @x[start[0]] + (val-bounds[0]) / df * dx
					ycross = @y[start[1]] + (val-bounds[0]) / df * dy
					ce[val][edge] = [xcross, ycross]
				end
				
			end
		end
		ce
	end
		
	private :get_crossed_edges
	
	def analyse_cell(i, j, cons, specific_val = nil)
		crossed_edges = get_crossed_edges(i, j, cons, specific_val)
		crossed_edges.each do |val, crossings|
			if specific_val
			 next	unless val == specific_val
			 #@analysed[val] ||= {}
			 #@analysed[val][[i,j]] = true
			 @found[[i,j]] ||= {}
			end
			outer = crossings.keys.find_all{|edge| edge < 4}
			inner = crossings.keys.find_all{|edge| edge > 4}
			next if outer.size == 0 and inner.size == 0
			VALID_CONNECTIONS.each do |connection|
				path = crossings.values_at(*connection).compact.zip(connection.map{|edge| [i, j, edge]})
# 					p path
				next if path.size != connection.size
				connection.each{|edge| crossings.delete(edge)}
				if specific_val
				 next if @found[[i,j]][path]
				 @found[[i,j]][path] = true

				end
				cons[val].push path
				
			end

# 				p val, crossings, cons
		end
	end
	
	private :analyse_cell
	
	def connect_contours(contours)
		return contours if contours.size == 1
		loop do
			catch(:restart) do 
# 				old_contours = contours
# 				contours = []
			joined = []
			for i in 0...contours.size
				break if i >= contours.size
				for j in i+1...contours.size
					break if j >= contours.size
					coni = contours[i]
					conj = contours[j]
					if joins?(coni[-1], conj[0])
						contours[i] = coni + conj
						contours.delete_at(j)
# 							redo
						throw(:restart)
					elsif joins?(coni[0], conj[-1])
						contours[i] = conj + coni
						contours.delete_at(j)
# 							redo
						throw(:restart)
					elsif joins?(coni[0], conj[0])
						contours[i] = coni.reverse + conj
						contours.delete_at(j)
# 							redo
						throw(:restart)
					elsif joins?(coni[-1], conj[-1])
						contours[i] = conj + coni.reverse
						contours.delete_at(j)
# 							redo
						throw(:restart)
					end


					
				end
			end
			contours.each{|con| con.uniq!}
			return contours
			end
			
		end
	end
	
	private :connect_contours
			
	
	
	def joins?(path, contour)	
# 			p @x.size
# 			dx = (@x.subvector(1, @x.size - 1) - @x.subvector(0, @x.size - 1)).sum / @x.size #The average dx between gridpoints
# 			dy = (@y.subvector(1, @y.size - 1) - @y.subvector(0, @y.size - 1)).sum / @y.size
# 			eps = 1.0
# # 			return true
# 			return ((path[0] - contour[0]).abs < dx.abs * eps and (path[1] - contour[1]).abs < dy.abs * eps)
# 			p dx; exit
		pi, pj, pedge = path[1]
		ci, cj, cedge = contour[1]
		joins = false
		joins = ((pi == ci and pj == cj + 2 and pedge == 0 and cedge == 2) or
			(pi == ci and pj == cj - 2 and pedge == 2 and cedge == 0) or
			(pi == ci + 2 and pj == cj  and pedge == 1 and cedge == 3) or
			(pi == ci - 2 and pj == cj  and pedge == 3 and cedge == 1))
# 			unless joins
# 				p path[1], contour[1]
# 				STDIN.gets
# 			end
# 			if path[1] == [2,10,1]
# 				p path[1], contour[1], joins
# 				STDIN.gets
# 			end
		return joins
	end
	
	private :joins?
		
		

	
end

module MultiFit
	class MultidLM
		def self.alloc(*args)
			new(*args)
		end
		
		def initialize(yproc = nil, fproc, dfproc, ndata, ndims, nparams)
			@fproc = Proc.new do |x, t, y, sigma, f|
# 				gridpoints = (0...@ndims).to_a.map do |i| 
# 					@gridpoints.col(i)
# 				end
				fproc.call(x, *@gridpoints, y, sigma, f)
			end
			@dfproc = Proc.new do |x, t, y, sigma, jac|
# 				gridpoints = (0...@ndims).to_a.map do |i| 
# 					@gridpoints.col(i)
# 				end
# 				puts 'hello'
				dfproc.call(x, *@gridpoints, y, sigma, jac)
			end

			@yproc = yproc
# 			fproc
# 			@dfproc = dfproc
			@ndata = ndata; @ndims = ndims; @nparams = nparams
			@f = GSL::MultiFit::Function_fdf.alloc(@fproc, @dfproc, @nparams)
			@solver = GSL::MultiFit::FdfSolver.alloc(FdfSolver::LMDER, @ndata, @nparams)
		end
		
		def set_data(xstart, *gridpoints, y, sigma)
# 			p 'g', gridpoints.size
			@gridpoints = gridpoints; @y = y; @x = xstart.dup; @sigma = sigma
			@t = GSL::Vector.alloc(@y.size)
			@t.set_all(0.0) # t should never be used.
			@f.set_data(@t, @y, @sigma)
			@solver.set(@f, @x)
		end
			
		def solve(print_out = false)
			(puts "Warning: due to a bug, print out doesn't work with less than 3 params"; print_out = false) if @nparams < 3
# 			p @nparams, @solver.send(:p)
			iter = 0
			@solver.print_state(iter) if print_out
			begin
			  iter += 1
			  status = @solver.iterate
			  @solver.print_state(iter) if print_out
			  status = @solver.test_delta(1e-7, 1e-7)
			end while status == GSL::CONTINUE and iter < 500
			
			@covar = @solver.covar(0.0)
			@position = @solver.position
			@my_chi2 = 0.0 
# 			gp = @gridpoints.transpose
			for i in 0...@y.size
				@my_chi2 += (@y[i] - eval(*@gridpoints.map{|vec| vec[i]}))**2.0 / @sigma[i]**2.0
			end
			@chi2 = (@solver.f.dnrm2)**2
			@dof = @ndata - @nparams
			@solved = true
			@solver.position
		end
		
		attr_accessor :chi2, :my_chi2, :covar, :position, :dof
		
		def eval(*points)
			raise "yproc not set" unless @yproc
			@yproc.call(@solver.position, *points)
		end

		
	end
	class MultidLMNumDiff < MultidLM
		def self.alloc(*params)
			new(*params)
		end

		def initialize(yproc, ndata, ndims, nparams)
			fproc = Proc.new do |x, *gridpoints, y, sigma, f|
# 				gridpoints = (0...@ndims).to_a.map do |i| 
# 					@gridpoints.col(i)
# 				end
				for i in 0...ndata.size do
					f[i] = (@yproc.call(x, *gridpoints.map{|vec| vec[i]}) - y[i])/sigma[i]
				end
			end
			dfproc = Proc.new do |x, *gridpoints, y, sigma, jac|
					for j in 0...nparams do
						xj = x[j]
						xplus = x.dup
						xplus[j] = xj + @delt[j]
						xminus = x.dup
						xminus[j] = xj - @delt[j]

						for i in 0...ndata do
							gp = gridpoints.map{|vec| vec[i]}
							yplus = @yproc.call(xplus, *gp)
							yminus = @yproc.call(xminus, *gp)
							#p "delt", @delt, "sigma", @sigma, "jac", jac, "jac.shape", jac.shape, "result", (yplus - yminus)/2*@delt[j]/sigma[i]
							jac.set(i, j, (yplus - yminus)/2*@delt[j]/sigma[i])
						end
					end
			end	
			super(yproc, fproc, dfproc, ndata, ndims, nparams)


		end
		def set_data(xstart, delt, *gridpoints, y, sigma)
			@delt = (delt || GSL::Vector.alloc([1e-5]*xstart.size))
			super(xstart, *gridpoints, y, sigma)
		end
	end


end



module SpectralAnalysis
	
# A Lomb periodogram is a method of spectrally analysing a set of data which are not evenly spaced, and thus cannot be Fourier transformed. The Lomb periodogram is something akin to a probability distribution function for a given set of frequencies.
	
class Lomb
	class << self
		alias :alloc :new
	end
	
	# Create a new Lomb object. times and data should be GSL::Vectors.
	
	def initialize(times, data)
		@times = times; @data = data
		raise "Times #{times.size} and data #{data.size} do not have the same size" unless @times.size == @data.size
		@n = data.size
		@dmean = data.mean
		@dvar = data.variance
	end
	
	attr_accessor :frequencies, :periodogram
	
	# Calculate the Lomb periodogram. Without studying the Lomb analysis, it's best to leave the defaults alone. frequency_factor is how far above the estimated Nyquist frequency (calculated by dividing the net time interval by the number of data points) the spectrum should be calculated. 
	
	def calculate_periodogram(frequency_factor=2.0, oversampling=4.0, frequency_indexes=nil)
		@frequency_factor = frequency_factor
		@oversampling = oversampling
		@nout = (@n * 0.5 * frequency_factor * oversampling).to_i # (nout or @n * 0.5 * frequency_factor * 4.0).to_i
		t_window = @times.max - @times.min
		delta_f = 1.0 / t_window / oversampling # / 2.0 / Math::PI #(@nout / @n / 0.5 / frequency_factor)
# 			data_min, data_max = @data.minmax
		
		@frequencies = GSL::Vector.linspace(delta_f, delta_f*@nout, @nout)
# 		p @nout, delta_f, @frequencies
		if frequency_indexes 
			@frequencies = GSL::Vector.alloc(frequency_indexes.map{|i| @frequencies[i]})
		end
		@periodogram = @frequencies.collect do |freq|
			p_n(freq)
		end
# 		@frequencies = @frequencies / Math::PI / 2.0
		[@frequencies, @periodogram]
	end
	
	# Proportional to the probability that the given frequency was present in the data. Roughly akin to p(k) for a Fourier transform.
	
	def p_n(freq)
		omega = freq * Math::PI * 2.0
		twoomt = @times * 2.0 * omega
		tau = Math.atan(
			twoomt.sin.sum / twoomt.cos.sum
			)/ 2.0 / omega
		omttau = ((@times - tau) * omega)
		c = omttau.cos
		s = omttau.sin
		ddmean = @data - @dmean
		pn = 1 / 2.0 / @dvar * (
			(ddmean * c).sum ** 2.0 / c.square.sum +
			(ddmean * s).sum ** 2.0 / s.square.sum
		)
		pn
	end
	
	# Equal to 1.0 - the probability that the value of pn could have been generated by gaussian noise
	
	def confidence(pn, frequency_factor = @frequency_factor)
		(1.0 - Math.exp(-pn)) ** (@n  * frequency_factor)
	end
	
	# The probability that the value of pn could have been generated by gaussian noise.
	
	def pnull(pn, frequency_factor = @frequency_factor)
		1.0 - confidence(pn, frequency_factor)
	end
		
	# Find a 
		
	def p_n_from_confidence(confidence, frequency_factor = @frequency_factor)
		- Math.log(1.0 - confidence ** (1.0 / @n / frequency_factor))
	end
		           
	
	def graphkit
		CodeRunner::GraphKit.autocreate(x: {title: "Frequency", data: @frequencies}, y: {title: "P_N", data: @periodogram})
	end
			
				
		
end
end

class GaussianSmoothKernel < Vector
	KERNELS = {}
	
	def self.alloc(sigma, delt = 1.0)
		return KERNELS[[sigma,delt]] if KERNELS[[sigma,delt]]
		npix ||= (3.0*sigma / delt).floor
		kernel = super(2*npix + 1)
		for i in 0...kernel.size
			j = (i - npix) * delt
			kernel[i] = Math.exp(- j**2 / 2.0 / sigma**2) / ( 2.0 * Math::PI * sigma**2)
		end
		KERNELS[[sigma,delt]] = kernel / kernel.sum
	end
	
end


class Vector
	def rectangular_smooth
		smooth = dup
		for i in 1...(self.size-1)
			smooth[i] = (self[i-1] + self[i] + self[i+1]) / 3.0
		end
		smooth
	end
	

	
	def gaussian_smooth!(sigma)
		new = gaussian_smooth(sigma)
		for i in 0...size do i
			self[i] = new[i]
		end
		return nil
	end
	def gaussian_smooth(sigma)
		npix = (3.0*sigma).floor
		smooth = dup
		smooth.set_all(0.0)
		kernel = GaussianSmoothKernel.alloc(sigma)# gaussian_smooth_kernel(sigma)

# 		p kernel
		for i in 0...smooth.size
			range = [([i - npix, 0].max), ([i + npix, smooth.size - 1].min)]
			ke = kernel.subvector(range[0] - i  + npix, range[1] - range[0] + 1)
			ke = kernel / ke.sum
			for j in range[0]..range[1]
				smooth[i] += self[j] * ke[j - i + npix]
			end
		end
		smooth
	end

	def mean_no_outliers(nstd=1.0)
		av = mean
		std = sd
		self.dup.delete_if{|val| (val - mean).abs > nstd * std}.mean
	end
# 				next if i + j < 0 or i + j >= smooth.size 
end

class Matrix
	def gaussian_smooth(sigmai, sigmaj = nil)
		sigmaj ||= sigmai
		for i in 0...shape[0]
			set_row i,row(i).gaussian_smooth(sigmai)
		end
		for i in 0...shape[1]
			set_col i,col(i).gaussian_smooth(sigmaj)
		end
	end
end


end


if ARGV.include? "test_coderunner_gsl_tools" and $has_put_startup_message_for_code_runner
# 
   if true 
 	
      x = []; y=[]; z=[]
 	
 	n = 10
 	zmat = GSL::Matrix.alloc(n, n)
 	for i in 0...n
 		for j in 0...n
 # 			next if rand < 0.5
 			x.push i; y.push j; z.push (i-0.0)**2 - (j-0.0)**3
 # 			x.push i; y.push j; z.push (i-2.0)**2 + (j-2.0)**2
 # 			x.push i; y.push j; z.push Math.exp((i-2.0)**2 + (j-2.0)**2 + 0.1)
 # 			x.push i * 18.0/n; y.push j* 18.0/n; z.push Math.sin((i* 18.0/n-0.0)**2 + (j* 18.0/n-0.0)**2)**2
 # 			x.push i * 18.0/n; y.push j* 18.0/n; z.push Math.sin(i* 18.0/n-0.0)**2
 
 			zmat[i,j] = z[-1]
 
 		end
 	end
 # 	p x, y, z
 	
 	xvec = GSL::Vector.alloc(x.uniq.sort)
 	yvec = GSL::Vector.alloc(y.uniq.sort)
 	
 	int = GSL::ScatterInterp.alloc(:thin_plate_splines, [x, y, z], 1.0)
 	p int.eval(2, 3), 4-27
 	p int.eval(8,7), 8**2-7**3
 	con = int.to_contour(60)
 # 	con = GSL::Contour.alloc(xvec, yvec, zmat)
 	contours = con.contours(10) #(-50, -40, -30, -20, -10, 0.0, 10)
 # 	p contour
 	graphs = contours.map do |val, cons|
 		unless cons[0]
 			nil
 		else
 			(cons.map do |con|
 				contour = con.transpose
 				kit = CodeRunner::GraphKit.autocreate({x: {data: contour[0]}, y: {data: contour[1], title: val.to_s}})
 				kit.data[0].with = "lp"
 				kit
 			end).sum
 		end
 	end
 	graphs.compact.reverse.sum.gnuplot({key: "off"})
 	end
 	
 	times = GSL::Vector.alloc(400.times.inject([0]){|a,i| a[i+1] = a[i] + rand/1.0; a })
 	data = (0.125*2.0*Math::PI*times).cos + (0.25*2.0*Math::PI*times).cos + times / times.max * 4.0
 	kit = CodeRunner::GraphKit.autocreate(x: {data: times}, y: {data: data})
 	kit.data[0].with = 'lp'
 	kit.gnuplot
 	lomb = GSL::SpectralAnalysis::Lomb.new(times, data)
 	lomb.calculate_periodogram(0.4)
 	kit = lomb.graphkit
 	kit.data[0].with = 'lp'
 	kit.gnuplot
	
	mat = GSL::Matrix.alloc(3, 3)
	mat.set_all 0.0
	mat[1,1] = 4.0
	mat.gaussian_smooth(0.8)
	p mat
	
	vec = GSL::Vector.linspace(0, 40, 100)
	vec = vec.gaussian_smooth(10.0)
	CodeRunner::GraphKit.quick_create([vec]).gnuplot
	
# 	STDIN.gets
	
	
end	
			
