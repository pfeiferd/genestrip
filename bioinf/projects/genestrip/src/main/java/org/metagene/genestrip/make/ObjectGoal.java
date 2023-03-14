package org.metagene.genestrip.make;

public abstract class ObjectGoal<T> extends Goal {
	private T object;
	
	public ObjectGoal(String name, Goal... dependencies) {
		super(name, dependencies);
	}
	
	public T get() {
		if (!isMade()) {
			make();
		}
		return object;
	}
	
	protected void set(T object) {
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Setting " + this + " to " + object);
		}
		this.object = object;
	}
	
	@Override
	public void cleanThis() {
		object = null;
	}

	@Override
	public boolean isMade() {
		return object != null;
	}
	
	@Override
	public String toString() {
		return "object goal: "+ getName();
	}
}
