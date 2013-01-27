/*!
	flux.js by morriswmz
	---------------------------------------
	A simple library to play with particles.
	Not very fancy, but works like a charm :)
	I would be appreciated if you give credit
	to me when you use it or extend it.
	
	NOTE : NO render is included
*/
(function () {
	var flux = {};
	// utils

	function extend(src, obj, obj2) {
		if (obj2) {
			// merge 3
			for (var key1 in obj2) {
				obj[key1] = obj2[key1];
			}
			for (var key2 in obj) {
				src[key2] = obj[key2];
			}
		} else {
			// merge 2
			for (var key in obj) {
				src[key] = obj[key];
			}
			return src;
		}
	}
	// add methods
	function methods(obj, fns) {
		if (typeof(obj) == 'function' ) {
			extend(obj.prototype, fns);
		} else {
			extend(obj, fns);
		}
	}
	// inherit from parent
	function inherit(child, parent) {
		var cons = child.prototype.constructor;
		if (typeof(parent.constructor) == 'function') {
			child.prototype = new parent();
			child.prototype._super = parent;
			child.prototype.constructor = cons;
		} else {
			child.prototype = parent;
			child.prototype._super = parent;
			child.prototype.constructor = cons;
		}
	}

	// geom & math
	flux.geom = {};
	// 2D Point
	flux.geom.Point = (function () {
		function Point(x, y) {
			this.x = x || 0.0;
			this.y = y || 0.0;
		}
		// static methods
		Point.distance = function (p1, p2) {
			var dx = p1.x - p2.x,
				dy = p1.y - p2.y;
			return Math.sqrt(dx*dx+dy*dy);
		}

		Point.Interpolate = function (p1, p2, t) {
			var dx = p2.x - p1.x,
				dy = p2.y - p1.y,
				nx = p1.x + dx * t,
				ny = p1.y + dy * t;
			return new Point(nx, ny);
		}
		// methods
		methods(Point, {
			clone : function () {
				return new Point(this.x, this.y);
			},
			equals : function (p) {
				return (this.x == p.x && this.y == p.y);
			},
			length : function () {
				return Math.sqrt(this.x * this.x + this.y * this.y);
			},
			normalize : function () {
				var l = this.length();
				if (l > 0) {
					this.x /= l;
					this.y /= l;
				}
			}
		});

		return Point;
	})();

	// Rectangle Class
	flux.geom.Rectangle = (function () {
		function Rectangle(x, y, width, height) {
			this.x = x;
			this.y = y;
			this.width = width;
			this.height = height;
		}

		methods(Rectangle, {
			top : function () { return this.y; },
			right : function () { return this.x + this.width; },
			bottom : function () { return this.y + this.height; },
			left : function () { return this.x; },
			contains : function (px, py) {
				return (px >= this.x &&
						px <= this.x + this.width &&
						py >= this.y &&
						py <= this.y + this.height);
			}
		});

		return Rectangle
	})();

	// particle definition
	flux.Particle = function () {
		this.x = 0.0;
		this.y = 0.0;
		this.z = 0.0;
		this.vx = 0.0;
		this.vy = 0.0;
		this.vz = 0.0;
		this.ax = 0.0;
		this.ay = 0.0;
		this.az = 0.0;

		this.lifespan = Number.MAX_VALUE;
		this.age = 0.0;

		this.properties = 0;
		// point to next particle
		this.next = null;
	};

	// particle manager
	flux.ParticleManager = (function (Particle, Rectangle) {
		// ts - timeStep 
		function PM(ts) {
			this.timeStep = ts || 1.0/30.0;
			this.activePool = new Particle();
			this.garbagePool = new Particle();
			this.activeCount = 0;
			this.garbageCount = 0;
			this.fields = [];
			this.boundary = new Rectangle(-32768.0, -32768.0, 65536.0, 65536.0);
		}
		// private functions
		function createParticles(n, lifespan) {
			var i = 0, list, tmp,
				self = this,
				gPool = self.garbagePool;
			// add head
			if (gPool.next) {
				list = gPool.next;
				gPool.next = list.next;
				list.next = null;
				list.age = 0;
				list.ax = 0;
				list.ay = 0;
				list.az = 0;
				list.lifespan = lifespan;
				self.garbageCount--;
			} else {
				list = new Particle();
				list.age = 0;
				list.lifespan = lifespan;
			}
			for (i = 1; i < n; i++) {
				if (gPool.next) {
					tmp = gPool.next;
					gPool.next = tmp.next;
					tmp.next = list.next;
					list.next = tmp;
					tmp.age = 0;
					tmp.ax = 0;
					tmp.ay = 0;
					tmp.az = 0;
					tmp.lifespan = lifespan;
					self.garbageCount--;
				} else {
					tmp = new Particle();
					tmp.next = list.next;
					list.next = tmp;
					tmp.age = 0;
					tmp.ax = 0;
					tmp.ay = 0;
					tmp.az = 0;
					tmp.lifespan = lifespan;
				}
			}
			self.activeCount += n;
			return list;
		}

		// methods
		methods(PM, {
			isActive : function () {
				return (this.activePool.next != null);
			},
			particleList : function () {
				return this.activePool.next;
			},
			update : function (fn) {
				var self = this,
					timeStep = self.timeStep,
					aPool = self.activePool,
					gPool = self.garbagePool,
					curParticle = aPool.next,
					prevParticle = aPool,
					tmp, i;
				while (curParticle) {
					curParticle.age += timeStep;
					if (curParticle.age >= curParticle.lifespan || !self.boundary.contains(curParticle.x, curParticle.y)) {
						// collect non-active particles
						tmp = curParticle.next;
						prevParticle.next = curParticle.next;
						curParticle.next = gPool.next;
						gPool.next = curParticle;
						curParticle.age = curParticle.lifespan;
						curParticle = tmp;
						self.garbageCount++;
						self.activeCount--;
					} else {
						curParticle.x += (curParticle.vx + 0.5 * curParticle.ax * timeStep) * timeStep;
						curParticle.y += (curParticle.vy + 0.5 * curParticle.ay * timeStep) * timeStep;
						curParticle.z += (curParticle.vz + 0.5 * curParticle.az * timeStep) * timeStep;
						curParticle.vx += curParticle.ax * timeStep;
						curParticle.vy += curParticle.ay * timeStep;
						curParticle.vz += curParticle.az * timeStep;
						// reset acceleration
						curParticle.ax = 0;
						curParticle.ay = 0;
						curParticle.az = 0;
					}
					if (curParticle) {
						prevParticle = curParticle;
						curParticle = curParticle.next;
					}
				}
				for (i = 0; i < self.fields.length; i++) {
					self.fields[i].apply(aPool.next);
				}
				if (fn != null) { fn.call(this, activePool.next) };
			},
			generate : function (emitter, n, lifespan, properties) {
				if (!emitter) throw new Error('No emitter specified');
				n = n || 0;
				lifespan = lifespan || Number.MAX_VALUE;
				properties = properties || 0;

				var newParticles = createParticles.call(this, n, lifespan),
					tail = emitter.initParticles(newParticles);
				tail.next = this.activePool.next;
				this.activePool.next = newParticles;
			},
			clear : function () {
				var self = this,
					aPool = self.activePool,
					gPool = self.garbagePool,
					curParticle = aPool.next,
					tmp;
				while (curParticle) {
					aPool.next = curParticle.next;
					tmp = curParticle.next;
					curParticle.next = gPool.next;
					gPool.next = curParticle;
					curParticle = tmp;
					self.garbageCount++;
					self.activeCount--;
				}
			}
		});

		return PM;
	})(flux.Particle, flux.geom.Rectangle);

	// Emitters
	flux.emitters = {};
	flux.emitters.EmissionDirection = {
		INSIDE : 0,
		OUTSIDE : 1,
		BOTH_SIDE : 2
	};
	flux.emitters.AbstractEmitter = function (v, vv) {
		this.velocityBase = v || 0.0;
		this.velocityVariance = vv || 0.0;
	}
	// Point Emitter
	flux.emitters.PointEmitter = (function () {
		function PointEmitter(config) {
			var self = this;
			self.centerX = config.centerX || 0.0;
			self.centerY = config.centerY || 0.0;
			self._super(config.velocityBase, config.velocityVariance);
			if (config.angleStart != undefined && config.angleStop != undefined) {
				self.angleStart = config.angleStart;
				self.angleStop = config.angleStop;
				if (self.angleStart < self.angleStop) {
					var tmp = self.angleStart;
					self.angleStart = self.angleStop;
					self.angleStop = tmp;
				}
			} else {
				self.angleStart = 0.0;
				self.angleStop = Math.PI * 2.0;
			}
		}

		inherit(PointEmitter, flux.emitters.AbstractEmitter);

		methods(PointEmitter, {
			setCenter : function (x, y) {
				this.centerX = x;
				this.centerY = y;
			},
			initParticles : function (particles) {
				var self = this,
					curParticle = particles,
					prevParticle,
					curV, curAngle,
					angleSpan = self.angleStop - self.angleStart;
				while (curParticle) {
					curParticle.x = self.centerX;
					curParticle.y = self.centerY;
					curV = self.velocityBase + (Math.random() - 0.5) * self.velocityVariance;
					curAngle = self.angleStart + Math.random() * angleSpan;
					curParticle.vx = curV * Math.cos(curAngle);
					curParticle.vy = curV * Math.sin(curAngle);
					prevParticle = curParticle;
					curParticle = curParticle.next;
				}
				return prevParticle;
			}
		});

		return PointEmitter;
	})();
	// Circular Emitter
	flux.emitters.CircularEmitter = (function (EmissionDirection) {
		function CircularEmitter(config) 
		{
			var self = this;
			self.centerX = config.centerX || 0.0;
			self.centerY = config.centerY || 0.0;
			self.radius = config.radius || 0.0;
			if (self.radius < 0.0) self.radius = 0.0;
			self._super(config.velocityBase, config.velocityVariance);
			if (config.angleStart != undefined && config.angleStop != undefined) {
				self.angleStart = config.angleStart;
				self.angleStop = config.angleStop;
				if (self.angleStart < self.angleStop) {
					var tmp = self.angleStart;
					self.angleStart = self.angleStop;
					self.angleStop = tmp;
				}
			} else {
				self.angleStart = 0.0;
				self.angleStop = Math.PI * 2.0;
			}
			self.direction = config.direction || EmissionDirection.OUTSIDE;
		}

		inherit(CircularEmitter, flux.emitters.AbstractEmitter);

		methods(CircularEmitter, {
			setCircle : function (x, y, r) {
				this.centerX = x;
				this.centerY = y;
				this.radius = r;
			},
			initParticles : function (particles) {
				var self = this,
					curParticle = particles,
					prevParticle,
					curV, curAngle,
					angleSpan = self.angleStop - self.angleStart,
					sinT, cosT,
					dirFactor = (self.direction == EmissionDirection.INSIDE) ? -1.0 : 1.0;
				switch (self.direction) {
					case EmissionDirection.BOTH_SIDE:
						while (curParticle) {
							dirFactor = Math.random() > 0.5 ? 1.0 : -1.0;
							curAngle = self.angleStart + Math.random() * angleSpan;
							sinT = Math.sin(curAngle);
							cosT = Math.cos(curAngle);
							curParticle.x = self.centerX + self.radius * cosT;
							curParticle.y = self.centerY + self.radius * sinT;
							curV = self.velocityBase + (Math.random() - 0.5) * self.velocityVariance;
							curParticle.vx = curV * cosT * dirFactor;
							curParticle.vy = curV * sinT * dirFactor;
							prevParticle = curParticle;
							curParticle = curParticle.next;
						}
						break;
					default:
						while (curParticle) {
							curAngle = self.angleStart + Math.random() * angleSpan;
							sinT = Math.sin(curAngle);
							cosT = Math.cos(curAngle);
							curParticle.x = self.centerX + self.radius * cosT;
							curParticle.y = self.centerY + self.radius * sinT;
							curV = self.velocityBase + (Math.random() - 0.5) * self.velocityVariance;
							curParticle.vx = curV * cosT * dirFactor;
							curParticle.vy = curV * sinT * dirFactor;
							prevParticle = curParticle;
							curParticle = curParticle.next;
						}
				}
				return prevParticle;
			}
		});
		
		return CircularEmitter;
	})(flux.emitters.EmissionDirection);

	// Path Emitter
	flux.emitters.PathEmitter = (function (EmissionDirection, Point) {		
		function PathEmitter(config) 
		{
			this._super(config.velocityBase, config.velocityVariance);
			dir = config.direction || EmissionDirection.BOTH_SIDE;
			this.setPathFromArray(config.path);
		}
		
		inherit(PathEmitter, flux.emitters.AbstractEmitter);

		methods(PathEmitter, {
			setPathFromArray : function (arr) {
				var self = this,
					i, n = arr.length;
				self._path = [];
				var path = self._path;
				// get points
				for (i = 0; i < n; i += 2) {
					path.push(new Point(arr[i], arr[i + 1]));
				}
				// update lengths
				self._cumulativeNormalizedLengths = [];
				self._normalizedLengths = [];
				self._lengths = [];
				var sum, // aliases
					cumulativeNormalizedLengths = self._cumulativeNormalizedLengths,
					normalizedLengths = self._normalizedLengths,
					lengths = self._lengths;
				n = path.length;
				sum = 0.0;
				for (i = 0; i < n - 1; i++) {
					normalizedLengths[i] = Point.distance(path[i], path[i + 1]);
					lengths[i] = normalizedLengths[i];
					sum += normalizedLengths[i];
					cumulativeNormalizedLengths.push(sum);
				}
				n = cumulativeNormalizedLengths.length;
				for (i = 0; i < n; i++) {
					cumulativeNormalizedLengths[i] /= sum;
					normalizedLengths[i] /= sum;
				}
			},
			initParticles : function (particles) {
				var self = this,
					curParticle = particles,
					prevParticle,
					curV, rnd, idx, n,
					sx, sy, ex, ey, dx, dy, l;
				while (curParticle) {
					// get segment index
					idx = 0;
					n = self._cumulativeNormalizedLengths.length
					rnd = Math.random();
					while (idx < n && rnd > self._cumulativeNormalizedLengths[idx]) {
						idx++;
					}
					if (idx > 0) rnd -= self._cumulativeNormalizedLengths[idx - 1];
					l = self._lengths[idx];
					rnd /= self._normalizedLengths[idx];
					// interpolate
					sx = self._path[idx].x;
					sy = self._path[idx].y;
					ex = self._path[idx + 1].x;
					ey = self._path[idx + 1].y;
					dx = ex - sx;
					dy = ey - sy;
					curParticle.x = sx + rnd * dx;
					curParticle.y = sy + rnd * dy;
					// update speed
					curV = (self.velocityBase + (Math.random() - 0.5) * self.velocityVariance);
					curV = (Math.random() > 0.5) ? curV : -curV;
					curParticle.vx = - curV * dy / l;
					curParticle.vy = curV * dx / l;
					prevParticle = curParticle;
					curParticle = curParticle.next;
				}
				return prevParticle;
			}
		});

		return PathEmitter;
	})(flux.emitters.EmissionDirection, flux.geom.Point);

	// fields
	flux.fields = {};
	// gravitation
	flux.fields.GravitationField = (function () {
		function GravitationField(gravityX, gravityY, gravityZ) 
		{
			this.gx = gravityX || 0.0;
			this.gy = gravityY || 0.0;
			this.gz = gravityZ || 0.0;
			
		}
		
		methods(GravitationField, {
			apply : function (particle) {
				var curParticle = particle;
				while (curParticle) {
					curParticle.ax += this.gx;
					curParticle.ay += this.gy;
					curParticle = curParticle.next;
				}
			}
		});

		return GravitationField;
	})();
	// random walker
	flux.fields.RandomWalkerField = (function () {
		function RandomWalkerField(angleRange) {
			this.angleRange = angleRange || Math.PI/12.0;
		}

		methods(RandomWalkerField, {
			apply : function(particle) {
				var curParticle = particle,
					dirVariance,
					cosT, sinT,
					tvx, tvy;
				while (curParticle) {
					dirVariance = this.angleRange * (Math.random() - 0.5);
					cosT = Math.cos(dirVariance);
					sinT = Math.sin(dirVariance);
					tvx = curParticle.vx;
					tvy = curParticle.vy;
					curParticle.vx = tvx * cosT - tvy * sinT;
					curParticle.vy = tvx * sinT + tvy * cosT;
					curParticle = curParticle.next;
				}
			}
		});

		return RandomWalkerField;
	})();
	// Resistance
	flux.fields.ResistanceField = (function () {
		function ResistanceField(resistanceFactor) 
		{
			this.k = resistanceFactor || 0.0;
			if (this.k < 0.0) this.k = 0.0;
		}
		
		methods(ResistanceField, {
			apply : function (particle) {
				var curParticle = particle;
				while (curParticle) {
					curParticle.ax += -this.k * curParticle.vx;
					curParticle.ay += -this.k * curParticle.vy;
					curParticle = curParticle.next;
				}
			}
		});

		return ResistanceField;
	})();
	// Repulsion
	flux.fields.RepulsionField = (function () {
		function RepulsionField(x, y, z, k) 
		{
			this.centerX = x || 0.0;
			this.centerY = y || 0.0;
			this.centerZ = z || 0.0;
			this._strength = k || 0.0;
			this._effectiveRadius = k;
		}
		
		methods(RepulsionField, {
			strength : function (val) {
				if (val) {
					if (val < 0) val = 0;
					this._strength = val;
					this._effectiveRadius = this._strength;
				} else {
					return this._strength;
				}
			},
			apply : function (particle) {
				var self = this,
					curParticle = particle,
					dx, dy, dist2, dist,
					fullForce;
				while (curParticle) {
					dx = curParticle.x - self.centerX;
					dy = curParticle.y - self.centerY;
					if (dx > self._effectiveRadius ||
						dx < -self._effectiveRadius ||
						dy > self._effectiveRadius ||
						dy < -self._effectiveRadius) {
						// out of effective range
						// no need to do complex calculation
						curParticle = curParticle.next;
					} else {
						dist2 = dx * dx + dy * dy + 1; // avoid divide by zero
						dist = Math.sqrt(dist2);
						fullForce = self._strength / dist2;
						curParticle.ax += fullForce * dx;
						curParticle.ay += fullForce * dy;
						curParticle = curParticle.next;
					}
				}
			}
		});
		
		return RepulsionField;
	})();
	// Attraction
	flux.fields.AttractionField = (function () {
		function AttractionField(x, y, z, k) 
		{
			this.centerX = x || 0.0;
			this.centerY = y || 0.0;
			this.centerZ = z || 0.0;
			this._strength = k || 0.0;
			this._effectiveRadius = Math.sqrt(k)*2.0;
		}
		
		methods(AttractionField, {
			strength : function (val) {
				if (val) {
					if (val < 0) val = 0;
					this._strength = val;
					this._effectiveRadius = Math.sqrt(val)*2.0;
				} else {
					return this._strength;
				}
			},
			apply : function (particle) {
				var self = this,
					curParticle = particle,
					dx, dy, dist2, dist,
					fullForce;
				while (curParticle) {
					dx = self.centerX - curParticle.x;
					dy = self.centerY - curParticle.y;
					if (dx > self._effectiveRadius ||
						dx < -self._effectiveRadius ||
						dy > self._effectiveRadius ||
						dy < -self._effectiveRadius) {
						// out of effective range
						// no need to do complex calculation
						curParticle = curParticle.next;
					} else {
						dist2 = dx * dx + dy * dy + 1; // avoid divide by zero
						dist = Math.sqrt(dist2);
						fullForce = self._strength / dist2 / dist;
						curParticle.ax += fullForce * dx;
						curParticle.ay += fullForce * dy;
						curParticle = curParticle.next;
					}
				}
			}
		});
		
		return AttractionField;
	})();

	// export to global namespace
	window.flux = flux;
})();