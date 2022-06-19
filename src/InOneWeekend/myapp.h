#pragma once
namespace Tmpl8
{

	class MyApp : public TheApp
	{
	public:
		// game flow methods
		void Init() {};
		void HandleInput() {};
		void Tick(float deltaTime) {};
		void Shutdown() { /* implement if you want to do something on exit */ }
		// input handling
		void MouseUp(int button) { mouseDown = false; }
		void MouseDown(int button) { mouseDown = true; }
		void MouseMove(int x, int y) { mousePos.x = x, mousePos.y = y; }
		void MouseWheel(float y) {};
		void KeyUp(int key) { /* implement if you want to handle keys */ }
		void KeyDown(int key) { /* implement if you want to handle keys */ }
		// data members
		float zoom = 100;							// map zoom
		int2 mousePos, dragStart, focusStart;		// mouse / map interaction
		bool mouseDown = false;						// keeping track of mouse button status
		static inline Kernel* kernel;
		//Buffer posBuffer;
		static inline cl_mem posDirBuffer,/*, dirBuffer,*/ peaksBuffer, frameBuffer, frameChangeBuffer;
		static inline uint counter = 0;
		static const inline uint particleSize = 7500;

		static inline float4 particlePosDir[particleSize]/*, particleDir[2 * particleSize]*/;
		static inline uint particleFrameChange[particleSize], particleFrame[particleSize], ids[particleSize];
		static inline uint count = 0;
	};

} // namespace Tmpl8
